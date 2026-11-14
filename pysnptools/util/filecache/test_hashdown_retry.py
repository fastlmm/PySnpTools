"""Tests for Hashdown download retries and atomic writes, using a fake urlopen."""

import io
import os
import tempfile
import unittest
from unittest import mock
from urllib.error import HTTPError, URLError

from pysnptools.util.filecache import hashdown
from pysnptools.util.filecache.hashdown import Hashdown


class _Body(io.BytesIO):
    """A response body; optionally fails partway through reading."""

    def __init__(self, data, fail_after=None):
        super().__init__(data)
        self._fail_after = fail_after

    def read(self, size=-1):
        if self._fail_after is not None and self.tell() >= self._fail_after:
            raise ConnectionResetError("simulated reset mid-download")
        return super().read(size if self._fail_after is None else min(size, 4))


def _http_error(code):
    return HTTPError("https://example.invalid/f", code, "simulated", {}, None)


class TestHashdownRetry(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.file = os.path.join(self.dir, "sub", "f.txt")
        sleep = mock.patch.object(hashdown.time, "sleep")
        self.sleep = sleep.start()
        self.addCleanup(sleep.stop)

    def _get(self, *responses):
        calls = []

        def fake_urlopen(url, timeout=None):
            calls.append(timeout)
            result = responses[len(calls) - 1]
            if isinstance(result, BaseException):
                raise result
            return result

        with mock.patch.object(hashdown, "urlopen", fake_urlopen):
            result = Hashdown._get_large_file("https://example.invalid/f", self.file, False)
        return result, calls

    def _leftovers(self):
        return [n for n in os.listdir(os.path.dirname(self.file)) if n.endswith(".part")]

    def test_retries_transient_network_error(self):
        result, calls = self._get(URLError("simulated DNS failure"), _Body(b"hello"))
        self.assertTrue(result)
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0], Hashdown._DOWNLOAD_TIMEOUT)
        with open(self.file, "rb") as f:
            self.assertEqual(f.read(), b"hello")
        self.sleep.assert_called_once_with(1)

    def test_retries_server_error(self):
        result, calls = self._get(_http_error(503), _http_error(502), _Body(b"ok"))
        self.assertTrue(result)
        self.assertEqual(len(calls), 3)
        self.assertEqual([c.args[0] for c in self.sleep.call_args_list], [1, 2])

    def test_404_returns_false_without_retry(self):
        result, calls = self._get(_http_error(404))
        self.assertFalse(result)
        self.assertEqual(len(calls), 1)
        self.sleep.assert_not_called()

    def test_client_error_is_not_retried(self):
        with self.assertRaises(HTTPError):
            self._get(_http_error(403))
        self.sleep.assert_not_called()

    def test_gives_up_after_all_attempts(self):
        errors = [URLError("down")] * Hashdown._DOWNLOAD_ATTEMPTS
        with self.assertRaises(URLError):
            self._get(*errors)
        self.assertEqual(self.sleep.call_count, Hashdown._DOWNLOAD_ATTEMPTS - 1)

    def test_interrupted_download_leaves_no_partial_file(self):
        os.makedirs(os.path.dirname(self.file), exist_ok=True)
        with open(self.file, "wb") as f:
            f.write(b"old complete file")
        result, calls = self._get(_Body(b"0123456789", fail_after=4), _Body(b"new"))
        self.assertTrue(result)
        self.assertEqual(len(calls), 2)
        with open(self.file, "rb") as f:
            self.assertEqual(f.read(), b"new")
        self.assertEqual(self._leftovers(), [])

    def test_failed_download_keeps_previous_file(self):
        os.makedirs(os.path.dirname(self.file), exist_ok=True)
        with open(self.file, "wb") as f:
            f.write(b"old")
        errors = [_Body(b"0123456789", fail_after=4) for _ in range(Hashdown._DOWNLOAD_ATTEMPTS)]
        with self.assertRaises(ConnectionResetError):
            self._get(*errors)
        with open(self.file, "rb") as f:
            self.assertEqual(f.read(), b"old")
        self.assertEqual(self._leftovers(), [])


def getTestSuite():
    return unittest.TestLoader().loadTestsFromTestCase(TestHashdownRetry)


if __name__ == "__main__":
    unittest.main()

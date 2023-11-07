import logging
import pysnptools.test
import unittest


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)

    suites = unittest.TestSuite([pysnptools.test.getTestSuite()])
    suites.debug

    r = unittest.TextTestRunner(failfast=True)  # cmk should be False
    ret = r.run(suites)
    assert ret.wasSuccessful()

    logging.info("done with testing")

"""Module sys to set python path to root"""
import sys

sys.path.append(".")

# pylint: disable=wrong-import-position
import libs.pbrpc as pbrpc  # noqa: E402

# pylint: enable=wrong-import-position


def test_pbrpc():
    # TODO

    assert True

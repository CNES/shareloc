"""Module sys to set python path to root"""
import sys

sys.path.append(".")

# pylint: disable=wrong-import-position
import libs.pbhelloworld as HWmodule  # noqa: E402

# pylint: enable=wrong-import-position


def test_helloworld():
    hw_object = HWmodule.HW()

    assert "Hello world !" == hw_object.hellow_world()

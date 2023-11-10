"""Module sys to set python path to root"""
import sys

sys.path.append(".")

# pylint: disable=wrong-import-position
import libs.pbhelloworld as HWmodule  # noqa: E402

# pylint: enable=wrong-import-position

HW_object = HWmodule.HW()
print(HW_object.hellow_world())

import sys
sys.path.append(".")
import libs.pbhelloworld as HWmodule


def test_helloworld():
    HWobject = HWmodule.HW()

    assert "Hello world !" == HWobject.hellow_world()
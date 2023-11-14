"""Module sys to set python path to root"""
import sys

sys.path.append(".")

# pylint: disable=wrong-import-position
import libs.pbrpc as pbrpc  # noqa: E402

# pylint: enable=wrong-import-position


row = [1.0, 1.0, 1.0]
col = [1.0, 1.0, 1.0]
alt = 1.0

GM_object = pbrpc.GeoModelTemplate("geomdel_path")
lonlatalt1 = GM_object.direct_loc_h(row, col, alt, False)
lonlatalt1 = GM_object.inverse_loc(row, col, alt)
del GM_object
# lonlatalt2 = GM_object.direct_loc_dtm(row,col,"test")
# lonlatalt3 = GM_object.inverse_loc(row,col,alt,)
# print(lonlatalt1)
# print(lonlatalt2)
# print(lonlatalt3)


# class RPCOptim(pbrpc.RPC):
#     def hello_python(self):
#         return "hello python"

# rpc_optim = RPCOptim()
# arg = [[1.,1.,1.],[1.,1.,1.]]
# zero  = rpc_optim.direct_loc_h(arg)

# print(type(zero))
# print(zero)

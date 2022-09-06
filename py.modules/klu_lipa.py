#
from lipa.solvers import *
try:
    #import hack_mkl
    import env
    from lipa_ctx import LIPA_ctx as lctx
    from lipa_ctx import LIPA_Green_z as lctx_z
    import lipa_ctx.time_iterators as itimes 
    LIPA_solver_group_ctx=lctx.LIPA_solver_group_ctx
    LIPA_solver_ctx=lctx.LIPA_solver_ctx
    LIPA_solver_ctx_z=lctx_z.LIPA_solver_ctx_z
    to_list=lctx_z.to_list
except Exception:
    from lipa.parallel import LIPA_solver_st
    LIPA_solver_ctx=LIPA_solver_st
    LIPA_solver_group_ctx=LIPA_solver_st
    import warnings
    warnings.warn('standart mode only [ LIPA_solver_st ] ')
    pass

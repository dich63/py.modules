

#include "mid2Dmesh.h"


using namespace mid2Dmesh;	
typedef matrix3x3_t<double>* pm3x3_t;

extern "C" void __cdecl  laplace_2D_to_matrix_trimesh
     (void*  pmatrix_proc,uint32_t index_base,uint32_t ntrs,uint32_t *ptrs,double *pa,double *pc
	 ,int64_t nvxs,double *pvxs,uint16_t * pfbc
	 ,double* pma,double* pmc,double* pmv)
{

	   
	to_matrix_proc_t matrix_proc=mesh_2D_matrix_procs_t<>::get_matrix_proc(pmatrix_proc);
	mesh_2D_t<>::laplace_2D_to_matrix_trimesh(matrix_proc,index_base,ntrs,ptrs,pa,pc,nvxs,pvxs,pfbc,pm3x3_t(pma),pm3x3_t(pmc),pm3x3_t(pmv));	
}

extern "C" void   vxs2D_prjtransform_in_index(uint32_t nvxs,uint32_t* pindex,double* pm,double* px,double* py)
{
	projective2D_t<>::vxs2D_prjtransform_in_index(nvxs,pindex,pm3x3_t(pm),px,py);
}
extern "C" void mesh2D_in_rect_index(uint32_t index_base,uint32_t fcenter,uint32_t ntrs,uint32_t* ptrs,uint32_t*pmasks
						  								 ,double* pvxs,uint32_t nrgn,uint32_t* prect_indx,double* prects)
{
    rect2D_t<double>* pr=(rect2D_t<double>*)prects;
	rect_utils_t<double>::mesh2D_in_rect_index(index_base,fcenter,ntrs,ptrs,pmasks,pvxs,nrgn,prect_indx,pr);
}
extern "C"  void vxs2D_in_rect_index(uint32_t fcenter,uint32_t nvxs,double* pvxs,uint32_t*pmasks
						 ,uint32_t nrgn,uint32_t* prect_indx,double* prects){
	rect2D_t<double>* pr=(rect2D_t<double>*)prects;
	rect_utils_t<double>::vxs2D_in_rect_index(fcenter,nvxs,pvxs,pmasks,nrgn,prect_indx,pr);
}



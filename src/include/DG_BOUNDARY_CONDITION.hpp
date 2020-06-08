#pragma once

struct DG_BOUNDARY_CONDITION
{
  double m_mem[16000];
  wmesh_int_t m_qn;
  wmesh_int_t m_nu;
  wmesh_int_t m_nx;
  wmesh_int_t m_ntest;

  
  //
  // Matrix to build the evaluation of u.
  //
  wmesh_mat_t<double> m_beval_uvw;
  wmesh_mat_t<double> m_beval_uvw_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  wmesh_mat_t<double> m_beval_xyz;
  wmesh_mat_t<double> m_beval_xyz_part[nfaceinelm];

  //
  // Matrix to build the evaluation of geometry of the element.
  //
  wmesh_mat_t<double> m_beval_test;
  wmesh_mat_t<double> m_beval_test_part[nfaceinelm];

  
  wmesh_mat_t<double> m_qrst;
  wmesh_mat_t<double> m_qrst_part[nfaceinelm];
  
  
  DG_BOUNDARY_CONDITION(){};
  DG_BOUNDARY_CONDITION(wmesh_shape_t* s_test_,
			wmesh_shape_t* s_teta_u_,
			wmesh_shape_t* s_teta_x_,
			wmesh_int_t 	qn,
			const double * __restrict__ 	qw,
			const double * __restrict__ 	qp)
  {
    const wmesh_int_t qnXnfaceinelm = qn * nfaceinelm;
    const wmesh_int_t nu = mkS_n(s_teta_u_);
    const wmesh_int_t nx = mkS_n(s_teta_x_);
    const wmesh_int_t ntest = mkS_n(s_test_);
    m_qn = qn;
    m_nu = nu;
    m_nx = nx;
    m_ntest = ntest;
    
    wmesh_int_t at = 0;

    //
    // Reference quadrature points 
    //
    {
      wmesh_int_t qrst_size = qnXnfaceinelm*dim;
      wmesh_mat_t<double>::define(&this->m_qrst,
			qnXnfaceinelm,
			dim,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&this->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    this->m_qrst.v + localFaceIndex * qn,
			    this->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   this->m_qrst.v,
		   &this->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      wmesh_mat_t<double>::define(&this->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&this->m_mem[at],
			nx);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&this->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    this->m_beval_xyz.v + (localFaceIndex * qn) * this->m_beval_xyz.ld,
			    this->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      wmesh_mat_t<double>::define(&this->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&this->m_mem[at],
			ntest);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&this->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    this->m_beval_test.v + (localFaceIndex * qn) * this->m_beval_test.ld,
			    this->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      wmesh_mat_t<double>::define(&this->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&this->m_mem[at],
			qnXnfaceinelm);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&this->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    this->m_beval_uvw.v + localFaceIndex * qn * this->m_beval_uvw.ld,
			    this->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    wmesh_int_t err;
    R rwork[4096];
    wmesh_int_t rwork_n = 4096;

    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      this->m_beval_xyz.v,
	      &this->m_beval_xyz.ld,
	      this->m_qrst.v,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      this->m_beval_uvw.v,
	      &this->m_beval_uvw.ld,
	      this->m_qrst.v,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      this->m_beval_test.v,
	      &this->m_beval_test.ld,
	      this->m_qrst.v,
	      &this->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  

    for (wmesh_int_t l=0;l<nfaceinelm;++l)
      {
	for (wmesh_int_t j=0;j<qn;++j)
	  {
	    for (wmesh_int_t i=0;i<ntest;++i)
	      {
		m_beval_test_part[l].v[m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }



  static void define(DG_BOUNDARY_CONDITION*self,
		     wmesh_shape_t* s_test_,
		     wmesh_shape_t* s_teta_u_,
		     wmesh_shape_t* s_teta_x_,
		     wmesh_int_t 	qn,
		     const double * __restrict__ 	qw,
		     const double * __restrict__ 	qp)
  {
    const wmesh_int_t qnXnfaceinelm = qn * nfaceinelm;
    const wmesh_int_t nu = mkS_n(s_teta_u_);
    const wmesh_int_t nx = mkS_n(s_teta_x_);
    const wmesh_int_t ntest = mkS_n(s_test_);
    self->m_qn = qn;
    self->m_nu = nu;
    self->m_nx = nx;
    self->m_ntest = ntest;
    
    wmesh_int_t at = 0;

    //
    // Reference quadrature points 
    //
    {
      wmesh_int_t qrst_size = qnXnfaceinelm*dim;
      wmesh_mat_t<double>::define(&self->m_qrst,
			qnXnfaceinelm,
			dim,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex < nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&self->m_qrst_part[localFaceIndex],
			    qn,
			    dim,
			    self->m_qrst.v + localFaceIndex * qn,
			    self->m_qrst.ld);
	}
      
      mkS_bmapping(qn,
		   self->m_qrst.v,
		   &self->m_qrst.ld,
		   qp);

      at += qrst_size;      
    }
    
#if 0
    matrix_handle_print(&m_qrst,stdout);
    printf("##\n");
#endif
    
    //
    // Evaluation of s_teta_x_
    //
    {      
      wmesh_mat_t<double>::define(&self->m_beval_xyz,
			nx,
			qnXnfaceinelm,
			&self->m_mem[at],
			nx);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&self->m_beval_xyz_part[localFaceIndex],
			    nx,
			    qn,
			    self->m_beval_xyz.v + (localFaceIndex * qn) * self->m_beval_xyz.ld,
			    self->m_beval_xyz.ld);
	}
      
      at += qnXnfaceinelm * nx;
    }

    //
    // Evaluation of s_test_
    //
    {      
      wmesh_mat_t<double>::define(&self->m_beval_test,
			ntest,
			qnXnfaceinelm,
			&self->m_mem[at],
			ntest);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&self->m_beval_test_part[localFaceIndex],
			    ntest,
			    qn,
			    self->m_beval_test.v + (localFaceIndex * qn) * self->m_beval_test.ld,
			    self->m_beval_test.ld);
	}
      
      at += qnXnfaceinelm * ntest;
    }

    //
    // Evaluation of s_teta_u_
    //    
    {      
      wmesh_mat_t<double>::define(&self->m_beval_uvw,
			nu,
			qnXnfaceinelm,
			&self->m_mem[at],
			qnXnfaceinelm);
      
      for (wmesh_int_t localFaceIndex=0;localFaceIndex<nfaceinelm;++localFaceIndex)
	{
	  wmesh_mat_t<double>::define(&self->m_beval_uvw_part[localFaceIndex],
			    nu,
			    qn,
			    self->m_beval_uvw.v + localFaceIndex * qn * self->m_beval_uvw.ld,
			    self->m_beval_uvw.ld);
	}
      
      at += qnXnfaceinelm * nu;
    }


    wmesh_int_t err;
    R rwork[4096*2];
    wmesh_int_t rwork_n = 4096*2;
    printf("-----------------------------------------------------\n");
    mkS_basis(mkS_b(s_teta_x_),
	      &qnXnfaceinelm,	      
	      self->m_beval_xyz.v,
	      &self->m_beval_xyz.ld,
	      self->m_qrst.v,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
    printf("-----------------------------------------------------\n");
    
    mkS_basis(mkS_b(s_teta_u_),
	      &qnXnfaceinelm,	      
	      self->m_beval_uvw.v,
	      &self->m_beval_uvw.ld,
	      self->m_qrst.v,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  
    printf("-----------------------------------------------------\n");

    
    mkS_basis(mkS_b(s_test_),
	      &qnXnfaceinelm,	      
	      self->m_beval_test.v,
	      &self->m_beval_test.ld,
	      self->m_qrst.v,
	      &self->m_qrst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);  
    printf("-----------------------------------------------------dddddddd\n");

    for (wmesh_int_t l=0;l<nfaceinelm;++l)
      {
	for (wmesh_int_t j=0;j<qn;++j)
	  {
	    for (wmesh_int_t i=0;i<ntest;++i)
	      {
		self->m_beval_test_part[l].v[self->m_beval_test_part[l].ld*j+i] *= qw[j];
	      }
	  }
      }
    
  }

struct DATA
{
  double m_mem[512];
  
  wmesh_mat_t<double> eval_xyz;
    wmesh_mat_t<double> eval_uvw;    

    vector_handle eval_udotn;
    vector_handle eval_f;
    vector_handle lrhs;

    DATA(wmesh_int_t qn,wmesh_int_t ntest)
    {
      // qn * (2 * dim + 3)
      wmesh_int_t at = 0;
      wmesh_mat_t<double>::define(&eval_xyz,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      wmesh_mat_t<double>::define(&eval_uvw,
			qn,
			dim,
			&m_mem[at],
			qn);
      at+=qn*dim;      
      vector_handle_def(&eval_udotn,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&eval_f,qn,&m_mem[at],1);
      at+=qn;
      vector_handle_def(&lrhs,ntest,&m_mem[at],1);
      at+=ntest;
    };
};    

DATA * CreateData()
{
  return new DATA(this->m_qn,
		  this->m_ntest);
};

  double m_eps;
  
  SmoothedHeaviside m_hea;
  template <typename userfct_t>
  void boundary_condition(const wmesh_int_t 		localFaceIndex,
			  const double * __restrict__ xu_,
			  const vector_handle&	normal,
			  const wmesh_mat_t<double>&	xyz,
			  const wmesh_mat_t<double>&	uvw,
			  DATA*			data,
			  userfct_t             userfct)
  {
    m_eps = 0.25;
    {
      Err e;
      SmoothedHeaviside_def	( &this->m_hea,
				  __eHeaviside_m4p3,
				  &m_eps,
				  &e);
    }
    
    
    //
    // Evaluate xyz.
    //
    data->eval_xyz = this->m_beval_xyz_part[localFaceIndex].transpose() * xyz;
#if 0    
    fprintf(stdout,"xyz " ifmt " " ifmt " " ifmt " " ifmt "\n",data->eval_xyz.n,data->eval_xyz.m,data->eval_xyz.ld,data->eval_f.ld);
    matrix_handle_print(&data->eval_xyz,stdout);
#endif

    //
    // User function.
    //


    //
    // eval_f pouvait ne pas etre initialise dans userfct
    //
    userfct(data->eval_xyz.n,
	    data->eval_xyz.v,
	    data->eval_xyz.ld,
	    data->eval_f.v,
	    data->eval_f.ld);
    
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif
    //
    // Evaluate uvw.
    //
#if 0
    fprintf(stdout,"uvw\n");
    matrix_handle_print(&uvw,stdout);
#endif
    
    data->eval_uvw = this->m_beval_uvw_part[localFaceIndex].transpose() * uvw;
    
#if 0
    fprintf(stdout,"eval_uvw\n");
    matrix_handle_print(&data->eval_uvw,stdout);
#endif    
    //
    // Compute the dot product with the normal
    //
    data->eval_udotn = data->eval_uvw * normal;
#if 0
    fprintf(stdout,"eval_udotn\n");
    vector_handle_print(&data->eval_udotn,stdout);
#endif

    //
    //
    //
    
    //
    // Form the flux functions.
    //
    
    for (wmesh_int_t i=0;i<data->eval_f.n;++i)
      {	//data->eval_f.v[i*data->eval_f.ld]=1.0;
	data->eval_f.v[i*data->eval_f.ld]
	  *= ( data->eval_udotn.v[i*data->eval_udotn.ld] < ((R)0.0) ) ? -data->eval_udotn.v[i*data->eval_udotn.ld] * xu_[0] : ((R)0.0) ;
      }
#if 0
    fprintf(stdout,"eval_f\n");
    vector_handle_print(&data->eval_f,stdout);
#endif	      
    //
    // Form the residu.
    //
    data->lrhs = this->m_beval_test_part[localFaceIndex] * data->eval_f;
#if 0
    fprintf(stdout,"lrhs\n");
    vector_handle_print(&data->lrhs,stdout);
#endif
  }

  

  void boundary_condition(ns_mesh * 	mesh,
			  const_wmesh_int_p  	cnc_u_,
			  const_wmesh_int_p  	cncoff_u_,
			  const double * __restrict__ 	data_u_,		       
			  const double * __restrict__	data_v_,
			  double * __restrict__ 		rhs_,
			  wmesh_int_t 		rhsoff_)
  {    

    vector_handle normal;
    wmesh_mat_t<double> xyz;
    wmesh_mat_t<double> uvw;
    double normal_values[6];
    double cooelm[32];
    double uvw_values[32];
    wmesh_mat_t<double>::define(&xyz,m_nx,dim,cooelm,m_nx);
    wmesh_mat_t<double>::define(&uvw,m_nu,dim,uvw_values,m_nu);
    vector_handle_def(&normal,2,normal_values,1);
    
    DATA* data = this->CreateData();
    for (wmesh_int_t jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {

		//
		//
		//
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const wmesh_int_t jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		normal_values[0]	= mesh->normaledge[2*jedge+0];
		normal_values[1] 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];
		//

		//
		//
		//
		
		for (wmesh_int_t k=0;k<m_nu;++k)
		  {
		    uvw_values[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];;
		  }
		for (wmesh_int_t k=0;k<m_nu;++k)
		  {
		    uvw_values[m_nu+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		  }
		//
		
		boundary_condition(jadj,
				   &xjac,
				   normal,
				   xyz,
				   uvw,
				   data,
				   [this](wmesh_int_t n_,const double * __restrict__ pos,wmesh_int_t posoff,double * __restrict__ f,wmesh_int_t foff)
				   {
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
#if 1
					 
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];
					 if (x<1.0e-13)
					   {
					     f[foff*k] = sin(12.0*y);
					   }
					 
					 
#else
					 
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];

					 
					 if ( (x < 1.0e-13) && (y>0.5 && y<=1.0) )
					   {
					     f[foff*k] = 0.0;
				       }
					 else if ( (x < 1.0e-13) && (y<=0.5 && y >0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x>0.5 && x<=1.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x<=0.5 && x >=0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (x > 1.0-1.0e-13) )
					   {
					 double z  = sin(acos(-1.0)*y);
					 f[foff*k] = z*z*z;
				       }
					 else
					   {
					 f[foff*k] = 0.0;
				       }
#endif
				       }
				     
				     
#if 0
				     
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = y - 0.5; // exp(-sin(32.0*y)*y);
				       }
				     Err e;
				     SmoothedHeaviside_eval	(&this->m_hea,
								 n_,
								 f,
								 &e);
#endif
#if 0
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = exp(-sin(32.0*y)*y);
				       }
#endif
				   });
		
		for (wmesh_int_t k=0;k<m_ntest;++k)
		  {
		    rhs_[jelm * rhsoff_ + k] += data->lrhs.v[k];
		  }		
	      }
	  }
      }
    
#if 0
    printf("ddddddddddddddddddd\n");
    for (wmesh_int_t jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (wmesh_int_t k=0;k<m_ntest;++k)
	  {
	    std::cout << rhs_[jelm * rhsoff_ + k] << std::endl;
	  }
	printf("\n");
      }
    printf("ddddddddddddddddddd\n");
#endif    

    
#if 0
    double mem[2048];
    wmesh_mat_t<double> beval_uvw[nfaceinelm];
    wmesh_mat_t<double> eval_uvw;
    wmesh_mat_t<double> beval_xyz[nfaceinelm];
    wmesh_mat_t<double> nrmelm;
    wmesh_mat_t<double> cooelm;
    wmesh_mat_t<double> uvwelm;
    vector_handle udotn;
    wmesh_mat_t<double> qxyz;
    vector_handle f;
  
    for (wmesh_int_t jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const wmesh_int_t jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		const R nx	= mesh->normaledge[2*jedge+0];
		const R ny 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];	    
		
	      //
	      // Interpolation of u * nx + v*ny + w * nz 
	      //

	      
	      //
	      // Evaluation of u on the edge.
	      //
	      eval_uvw = beval_uvw[jadj] * dofs_uvw;

	      //
	      // Evaluation of u on the edge.
	      //
	      udotn = eval_uvw * nrmedge;

	      //
	      //
	      // Cartesian coordinates of the quadrature points.
	      //
	      //
	      qxyz = build_qxyz * cooelm;

	      //
	      // Evaluation of the user function.
	      //
	      

	      //
	      // Form the flux functions.
	      //
	      for (wmesh_int_t i=0;i<nq;++i)
		{
		  f.v[i] *= ( udotn.v[i] < ((R)0.0) ) ? -udotn.v[i] * xu_[0] : ((R)0.0) ;
		}
	      
	      //
	      // Form the residu.
	      //
	      lrhs = qbasis * f;

	    }
	}
    }
#endif
  };

  
  void boundary_condition(CG_VAR& 	velocity_,
			  DG_VAR&       residual_)
  {    
    ns_mesh * mesh = velocity_.m_mesh;
    vector_handle normal;
    wmesh_mat_t<double> xyz;
    wmesh_mat_t<double> uvw;
    double normal_values[6];
    double cooelm[32];
    double uvw_values[32];
    wmesh_mat_t<double>::define(&xyz,m_nx,dim,cooelm,m_nx);
    wmesh_mat_t<double>::define(&uvw,m_nu,dim,uvw_values,m_nu);
    vector_handle_def(&normal,2,normal_values,1);
    
    DATA* data = this->CreateData();
    for (wmesh_int_t jelm=0;jelm<mesh->nelm;++jelm)
      {
	for (int jadj=0;jadj<3;++jadj)
	  {
	    if (mesh->adj[jelm*3+jadj] == 0)
	      {

		//
		//
		//
		ns_mesh_cooelm(mesh,
			       &jelm,
			       cooelm);
		
		const wmesh_int_t jedge   	= mesh->cnc[6*jelm+3+jadj]-mesh->m_numEntities[0]-1;  
		normal_values[0]	= mesh->normaledge[2*jedge+0];
		normal_values[1] 	= mesh->normaledge[2*jedge+1];
		const R xjac     	= mesh->jacedge[jedge];
		//

		//
		//
		//
		velocity_.dofselm(jelm, 0, uvw_values, 1);
		velocity_.dofselm(jelm, 1, &uvw_values[m_nu], 1);
#if 0
		for (wmesh_int_t k=0;k<m_nu;++k)
		  {
		    uvw_values[k] = data_u_[cnc_u_[jelm*cncoff_u_[0]+k]-1];;
		  }
		for (wmesh_int_t k=0;k<m_nu;++k)
		  {
		    uvw_values[m_nu+k] = data_v_[cnc_u_[jelm*cncoff_u_[0]+k]-1];
		  }
#endif
		//
		
		boundary_condition(jadj,
				   &xjac,
				   normal,
				   xyz,
				   uvw,
				   data,
				   [this](wmesh_int_t n_,const double * __restrict__ pos,wmesh_int_t posoff,double * __restrict__ f,wmesh_int_t foff)
				   {
				     // std::cout << "----" << std::endl;
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
#if 1
					 double x = pos[posoff*0+k];
					 double y = pos[posoff*1+k];
#if 0
					 if (x<1.0e-10)
					   {
					     f[foff*k] = sin(12.0*y);
					   }
#endif
#if 1
					 //					 std::cout << " " << x << " " << y << std::endl;
					 if ( (x < 1.0e-13) && (y>0.5 && y<=1.0) )
					   {
					 f[foff*k] = 0.0;
				       }
					 else if ( (x < 1.0e-13) && (y<=0.5 && y >0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x>0.5 && x<=1.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (y < 1.0e-13) && (x<=0.5 && x >=0.0) )
					   {
					 f[foff*k] = 1.0;
				       }
					 else if ( (x > 1.0-1.0e-13) )
					   {
					 double z  = sin(acos(-1.0)*y);
					 f[foff*k] = z*z*z;
				       }
					 else
					   {
					     f[foff*k] = 0.0;
				       }
#endif
#endif

				       }
#if 0
				     
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = y - 0.5; // exp(-sin(32.0*y)*y);
				       }
				     Err e;
				     SmoothedHeaviside_eval	(&this->m_hea,
								 n_,
								 f,
								 &e);
#endif
				     
#if 0
				     for (wmesh_int_t k=0;k<n_;++k)
				       {
					 double y = pos[posoff*1+k];
					 f[k] = exp(-sin(32.0*y)*y);
				       }
#endif
				   });

		//		residual_.setdofselm(jelm,0,data->lrhs.v,1);
		WLA::matrix_h r = residual_.matrix(); 
		for (wmesh_int_t k=0;k<m_ntest;++k)
		  {
		    //		    residual_.m_values3.v[jelm * residual_.m_values3.ld + k] += data->lrhs.v[k];
		    r.v[jelm * r.ld + k] += data->lrhs.v[k];
		  }		
	      }
	  }
      }
    
  }

};

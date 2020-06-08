#pragma once

struct DG_VIEW
{
  static constexpr I dim = 2;
  //  Generator<__eTopology_TRIANGLE,_degree> *generator;
  Uniform<__eTopology_TRIANGLE,5> *generator;


  I 		m_nsubcells;
  R 		m_mem[4096];

  
  WLA::matrix_h m_rst;
  WLA::matrix_h m_beval_xyz;
  WLA::matrix_h m_beval_f;
  WLA::matrix_h m_xyz;
  WLA::vector_h m_f;

  R p[32];

  DG_VIEW(mkS shape_x,mkS shape_f)
  {
    //Uniform<__eTopology_TRIANGLE,4> treilli;(p);
    p[0]=-2.0/3.0;
    p[1]=-1.0/3.0;
    p[2]=0;
    p[3]=1.0/3.0;
    p[4]=2.0/3.0;
#if 0
    p[0]=-0.9061798459386639927976268;
    p[1]=    -0.5384693101056830910363;
    p[2]=0;
    p[3]=0.5384693101056830910363144207;
    p[4]=0.9061798459386639927976268;
#endif

#if 0
p[0] =     -0.973906528517171720077964;
p[1] =     -0.8650633666889845;
p[2] =     -0.679409568299024;
p[3] =     -0.4333953941292471907992659;
p[4] =     -0.14887433898163121088482600112971;
p[5] =     0.1488743389816312108848260011;
p[6] =     0.433395394129247190799;
p[7] =     0.67940956829902440;
p[8] =     0.865063366688984;
p[9] =     0.9739065285171717200779640;
#endif


#if 0
    p[0] = -0.774596669241483377035853079956;
    p[1] = 0;
    p[2] = 0.774596669241483377035853079956;
    p[0] = -0.5;
    p[1] = 0;
    p[2] = 0.6;
#endif
    
    // generator = new Generator<__eTopology_TRIANGLE,_degree> (p);
    
    generator 	= new Uniform<__eTopology_TRIANGLE,5>();
    I nx 	= mkS_n(shape_x);
    I n 	= mkS_n(shape_f);
    I nsubcells = generator->nsubcells();    
    I at     	= 0;
    I nnodes 	= generator->nnodes();
    m_nsubcells = nsubcells;
    
    WLA::matrix_h::define(m_rst,nnodes,dim,&m_mem[at],nnodes);
    at += nnodes*dim;
    WLA::matrix_h::define(m_beval_xyz,nx,nnodes,&m_mem[at],nx);
    at += nnodes*nx;
    WLA::matrix_h::define(m_beval_f,n,nnodes,&m_mem[at],n);
    at += n*nnodes;
    WLA::matrix_h::define(m_xyz,nnodes,dim,&m_mem[at],nnodes);
    at += nnodes*dim;
    WLA::vector_h::define(m_f,nnodes,&m_mem[at],1);
    at += nnodes;

    for (I i=0;i<nnodes;++i)
      {
	m_rst.x[i] = generator->GetCoordinate<double>(i,0);
      }
    for (I i=0;i<nnodes;++i)
      {
	m_rst.x[m_rst.ld+i] = generator->GetCoordinate<double>(i,1);
      }

    I err;
    R rwork[4096];
    I rwork_n = 4096;
    
    mkS_basis(mkS_b(shape_x),
	      &nnodes,	      
	      this->m_beval_xyz.x,
	      &this->m_beval_xyz.ld,
	      this->m_rst.x,
	      &this->m_rst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);

    mkS_basis(mkS_b(shape_f),
	      &nnodes,	      
	      this->m_beval_f.x,
	      &this->m_beval_f.ld,
	      this->m_rst.x,
	      &this->m_rst.ld,	      
	      rwork,
	      &rwork_n,
	      &err);
  };
  
  void xyz(I subcell,
	   WLA::matrix_h&xyz)
  {
    for (I j=0;j<dim;++j)
      {
 	for (I i=0;i<3;++i)
	  {
	    xyz.x[xyz.ld*j+i] = m_xyz.x[ m_xyz.ld * j + generator->GetNodeIndex(subcell,i)];
	  }
      }
  };
  
  void f(I subcell,
	 WLA::vector_h&f)
  {
    for (I i=0;i<3;++i)
      {
	f.x[f.ld*i] = this->m_f.x[ this->m_f.ld * generator->GetNodeIndex(subcell,i)];
      }
    //    vector_handle_print(&f,stdout);
    //    exit(1);
  };
  
  void xyz(const WLA::matrix_h&xyz)
  {
    //    std::cout << "#############################################" << std::endl;
    //    std::cout << xyz << std::endl;
    //    std::cout << "#############################################" << std::endl;
    //    std::cout << this->m_beval_xyz << std::endl;
    //    std::cout << "#############################################" << std::endl;
    this->m_xyz = this->m_beval_xyz.transpose() * xyz;
    //  std::cout << "#############################################" << std::endl;
    
    //    std::cout << this->m_xyz << std::endl;
#if 0
    double r1=1;
    double r0=0;
    BlasMKL::gemm("T",
		  "N",
		  &m_xyz.m,
		  &m_xyz.n,
		  &xyz.n,
		  &r1,
		  this->m_beval_xyz.x,
		  &this->m_beval_xyz.ld,
		  xyz.x,
		  &xyz.ld,
		  &r0,
		  m_xyz.x,
		  &m_xyz.ld);

    std::cout << this->m_xyz << std::endl;
    exit(1);
#endif
    
  };

  void f(const WLA::vector_h&f)
  {
    this->m_f = this->m_beval_f.transpose() * f;
  };
  
};

#pragma once 
#include <string.h>  
#include <stdio.h>
namespace WCOMMON
{
  //!
  //!   @brief Class to define a command line.
  //!
  class cmdline
  {
  
  public:
    static constexpr int s_str_len = 256;
    typedef char str_t[s_str_len];
    
    //!
    //! @brief Define the command line
    //! @param argc_ nb arguments
    //! @param argv_ array of arguments
    //! 
    cmdline(int	argc_,
	    char** 	argv_);
  
    //!  
    //! @brief Free the cmdline object
    //!  
    ~cmdline();

    
    //!  
    //! @brief Display the command line
    //!   
    void disp() const;

    //!  
    //! @brief Display the command line
    //! @param out_ the file descriptor.
    //!   
    void disp_header(FILE*out_) const;
  
    //!  
    //! @brief Get a logical option
    //! @param opt_ name of the option
    //! @return false if not found, true otherwise
    //!  
    bool	get_logical(const char 	*opt_);

    //!  
    //! @brief Get a string option
    //! @param opt_ name of the option
    //! @param v_ value of the option
    //! @return false if not found, true otherwise
    //!  
    bool get_string	(const char* 	opt_,
			 char*   	v_);

    bool option(const char* 	opt_);
        
    template<typename T>
    bool option(const char* 	opt_,
		T*		v_);
    
    //!  
    //! @brief Get a real option
    //! @param opt_ name of the option
    //! @param v_ value of the option
    //! @return false if not found, true otherwise
    //!  
    template<typename real_t>
    bool  		get_real	(const char* 	opt_,
					 real_t*	v_);
    
    //!  
    //! @brief Get an integer option
    //! @param opt_ name of the option
    //! @param v_ value of the option
    //! @return false if not found, true otherwise
    //!
    template<typename int_t>
    bool  		get_integer	(const char* 	opt_,
					 int_t*		v_);
    
    //!  
    //! @brief Get the ith argument of the command line
    //! @param ith_  ith argument
    //! @return the ith argument 
    //!  
    const char* 		get_arg		(const int 	ith_) const;

    //!  
    //! @brief Get the number of arguments.
    //! @return The number of arguments.
    //!  
    inline int 			get_nargs	() const;

    //!  
    //! @brief Check invalid argument
    //! @return 1 if an invalid argument is found, 0 otherwise
    //!  
    bool 			check_invalid() const;
    
    //!  
    //! @brief Is empty ?
    //! @return false if not empty, true otherwise
    //!  
    bool			isempty() const;

  template<typename T>
  static void scan(const char* 	opt_,
	    T*			v_);

  private:
  
    static constexpr int MAXLEN = 256;
    static constexpr int MAXARG = 128;

    int 		argc{};			//!< number of arguments in command line 
    char 		argv[MAXARG][MAXLEN]; //!< arguments in command line 
    int 		ref_argc{};		//!< reference number of argument in command line 
    char **	ref_argv{};		//!< reference arguments in command line
  
    int search(int 		argc_,
	       char argv_[cmdline::MAXARG][cmdline::MAXLEN] ,
	       const char* 	opt_,
	       const int 	nchoice_ = 0,
	       const char** opt_choice_ = NULL);
  
  };

#include <string.h>  
#include <stdio.h>  



  //!  
  //! @brief Get the number of arguments.
  //! @return The number of arguments.
  //!  
  inline int cmdline::get_nargs	() const { return this->argc;}

  //!
  //!   @brief Class to define a command line.
  //!
  int cmdline::search(int 		argc_,
		      char argv_[cmdline::MAXARG][cmdline::MAXLEN] ,
		      const char* 	opt_,
		      const int 	nchoice_,
		      const char** opt_choice_)
  {
    int j,i = (int)0;
    for (i=1;i<argc_;i++)
      {
	if (!strcmp(argv_[i],opt_))
	  {
	    break;
	  }
      }
    if (i<argc_)
      {
	if ( (nchoice_>0) && (i+1<argc_) )
	  {
	    ++i;
	    for (j=0;j<nchoice_;++j)
	      {
		if (! strcmp(argv_[i],opt_choice_[j]))
		  {
		    break;
		  }
	      }
	    if (j<nchoice_)
	      {
		return i-1;
	      }
	    else 
	      {
		return (int)-1;
	      }
	  }      
	else
	  {
	    return i;
	  }
      }
    else
      {
	return (int)-1;
      }
  };

  void cmdline::disp_header(FILE * out_) const
  {

    fprintf(out_,"// date    : %s\n",__DATE__);
    fprintf(out_,"// command :");
    for (int i=0;i<this->argc;++i) 
      { 
	fprintf(out_," %s",this->argv[i]);
      } 
    fprintf(out_,"\n");
  };

  void cmdline::disp() const
  {

    for (int i=0;i<this->argc;++i) 
      { 
	fprintf(stdout,"[%d] %s\n",i,this->argv[i]);
      } 
  };
  
  cmdline::~cmdline()
  {
  };

  cmdline::cmdline(int		argc_,
		   char ** 		argv_)
  {  
    if (argc_<cmdline::MAXARG)
      {
	this->ref_argc 	= argc_;
	this->ref_argv 	= argv_;
	this->argc 	= argc_;
	for (int i=0;i<argc_;++i) 
	  { 
	    strncpy( this->argv[i],argv_[i],(size_t)cmdline::MAXLEN); 
	  } 
      }
    else
      {
	// ns_errmsg("cmdline::def:__mcmdline::MAXARG=%d is reached",argc_);
      }
  };



  bool cmdline::get_logical(const char*	opt_)
  {  
    int k = cmdline::search(this->argc,
			    this->argv,
			    opt_);
  
    if ( (k>=0) && (k<this->argc) )
      {
	if (k+1<this->argc)
	  {
	    { int i;
	      for (i=0;i<this->argc-k;++i)
		{
		  strncpy(this->argv[k+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
		} }
	  }
	this->argc-=1;
	return true;
      }
    return false;
  };

  bool cmdline::option(const char*	opt_)
  {  
    int k = cmdline::search(this->argc,
			    this->argv,
			    opt_);
  
    if ( (k>=0) && (k<this->argc) )
      {
	if (k+1<this->argc)
	  {
	    { int i;
	      for (i=0;i<this->argc-k;++i)
		{
		  strncpy(this->argv[k+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
		} }
	  }
	this->argc-=1;
	return true;
      }
    return false;
  };


  bool cmdline::get_string(const char* 		opt_,
			   char* 			t_)
  {  
    int j,k = cmdline::search(this->argc,
			      this->argv,
			      opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j=k-1;
	strncpy(t_,this->argv[k],(size_t)cmdline::MAXLEN);
	{
	  int i;
	  for (i=0;i<this->argc-k;++i)
	    {
	      strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	    }
	}
	this->argc-=2;
	return true;
      }
    return false;
  };



  template<typename int_t>
  bool  cmdline::get_integer(const char*	opt_,
			     int_t* 	x_)
  {  
    int j,k = cmdline::search(this->argc,
			      this->argv,
			      opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j=k-1;
	sscanf(this->argv[k],"%d",x_);
	{	int i;
	  for (i=0;i<this->argc-k;++i)
	    {
	      strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	    } }
	this->argc-=2;
	return true;
      }
    return false;
  };

  template<>
  bool  cmdline::get_integer<long long int>(const char*	opt_,
					    long long int* 	x_)
  {  
    int j,k = cmdline::search(this->argc,
			      this->argv,
			      opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j=k-1;
	sscanf(this->argv[k],"%Ld",x_);
	{	int i;
	  for (i=0;i<this->argc-k;++i)
	    {
	      strncpy(this->argv[j+i],this->argv[k+1+i],(size_t)cmdline::MAXLEN);
	    } }
	this->argc-=2;
	return true;
      }
    return false;
  };
  

  template<>
  void cmdline::scan(const char* 	opt_,
			    float*		v_)
  {
    sscanf(opt_, "%e", v_);
  }

  template<>
  void cmdline::scan(const char* 	opt_,
	    int*		v_)
  {
    sscanf(opt_, "%d", v_);
  }

  template<>
  void cmdline::scan(const char* 	opt_,
	    long int*		v_)
  {
    sscanf(opt_, "%ld", v_);
  }

  template<>
  void cmdline::scan(const char* 	opt_,
	    long long int*	v_)
  {
    sscanf(opt_, "%Ld", v_);
  }

  template<>
  void cmdline::scan(const char* 	opt_,
	    char*	v_)
  {
    sscanf(opt_, "%s", v_);
  }

  template<typename T>
  bool cmdline::option(const char* 	opt_,
		       T*		v_)
  {
    int j, k = cmdline::search(this->argc,
			       this->argv,
			       opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j = k-1;
	cmdline::scan(this->argv[k], v_);
	for (int i=0;i<this->argc-k;++i)
	  {
	    strncpy(this->argv[j+i],
		    this->argv[k+1+i],
		    size_t(cmdline::MAXLEN));
	  } 
	this->argc-=2;
	return true;
      }
    return false;
  }

  template<typename real_t>
  bool  cmdline::get_real(const char* 		opt_,
			  real_t* 		x_)
  {  
    int j,k = cmdline::search(this->argc,
			      this->argv,
			      opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j=k-1;
	sscanf(this->argv[k],"%le",x_);
	{	int i;
	  for (i=0;i<this->argc-k;++i)
	    {
	      strncpy(this->argv[j+i],this->argv[k+1+i],size_t(cmdline::MAXLEN));
	    } }
	this->argc-=2;
	return true;
      }
    return false;
  };


  template<>
  bool  cmdline::get_real<double>(const char* 		opt_,
				  double* 		x_)
  {  
    int j,k = cmdline::search(this->argc,
			      this->argv,
			      opt_);
    if ( (k>=0) && (++k<this->argc) )
      {
	j=k-1;
	sscanf(this->argv[k],"%le",x_);
	{	int i;
	  for (i=0;i<this->argc-k;++i)
	    {
	      strncpy(this->argv[j+i],this->argv[k+1+i],size_t(cmdline::MAXLEN));
	    } }
	this->argc-=2;
	return true;
      }
    return false;
  };


  bool cmdline::check_invalid() const
  {
    int i,j=0;
    for (i=2;i<this->argc;++i)
      {
	if (this->argv[i][0]=='-')
	  {
	    j=1;
	  }
      }
    return j>0; 
  };


  bool cmdline::isempty() const
  {
    return (1 == this->argc);
  };

  const char*	cmdline::get_arg(const int i_) const 
  {
    return (i_<this->argc) ? (this->argv[i_]) : nullptr;
  };
};

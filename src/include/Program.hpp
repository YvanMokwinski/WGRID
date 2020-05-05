#pragma once

//#include "System/Type.h"
//#include "Type.h"
#include "CmdLine.h"
#include "Mpi.h"
#include "Logging/Error.hpp"
#include "Logging/Warning.hpp"
#include <list>
#include <stdarg.h>
class OptionBase
{
protected:
  OptionBase(){};
public:
  virtual ~OptionBase(){};
  virtual void Analyze(CmdLine * cmdLine_) = 0;
  virtual void Usage() const = 0;
  virtual void CheckError()const = 0;
  virtual void CheckWarning()const = 0;
};

template <typename _type> struct TypeTypeString
{
  static std::string Get();
};

template <> struct TypeTypeString<int>
{
  static std::string Get()
  {
    return std::string("<integer>");
  };
};

template <> struct TypeTypeString<long int>
{
  static std::string Get()
  {
    return std::string("<integer>");
  };
};


template <> struct TypeTypeString<long long int>
{
  static std::string Get()
  {
    return std::string("<integer>");
  };
};

template <> struct TypeTypeString<double>
{
  static std::string Get()
  {
    return std::string("<real>");
  };
};

template <> struct TypeTypeString<std::string>
{
  static std::string Get()
  {
    return std::string("<string>");
  };
};

template <> struct TypeTypeString<bool>
{
  static std::string Get()
  {
    return std::string("");
  };
};

template <typename _type> class Option : public OptionBase
{
protected:
  std::string m_token;
  std::string m_description;
  _type m_defaultValue;
  _type *m_value;

  bool 	m_activate;
  bool m_isOptional;

public:

  virtual   std::string GetToken() const { return m_token; };
  virtual   bool IsActive() const { return m_activate; };
  virtual   bool IsOptional() const { return m_isOptional; };
  virtual   _type GetDefaultValue() const { return *m_value; };
  
 Option(_type * valuePointer_,
	const std::string& token_,
	const _type& defaultValue_,
	const bool& isOptional_,
	const std::string&description_ )

   : m_token(token_),
  m_description(description_),  
  m_defaultValue(defaultValue_),
  m_value(valuePointer_),
  m_activate(false),
  m_isOptional(isOptional_)
  {
    *m_value = m_defaultValue;
  };
 
 virtual void CheckError() const
 {
   if (!this->IsActive()) 
     {
       if (!this->IsOptional())
	 {
	   std::string msg = "{" + this->GetToken() + " " + TypeTypeString<_type>::Get() +"}" + std::string(" is not optional");
	   Logging::Error::Get() << msg << std::endl;
	   exit(1);
	   //	   throw(Exception(msg));
	 }
     }
 };

 virtual void CheckWarning() const
 {
   if (!this->IsActive()) 
     {
       Logging::Warning::Get() << "option '" << this->GetToken() <<"' set to default_value "<< m_defaultValue << std::endl; 
     }
 };

 virtual void Usage() const
 {
   if (m_isOptional)
     {
       std::cerr << "["  << m_token << " " << TypeTypeString<_type>::Get() <<"]" ;
     }
   else
     {
       std::cerr << "{"  << m_token << " " << TypeTypeString<_type>::Get() <<"}" ;
     }
   std::cerr <<  std::endl << "\t\t" << "default value: "<< m_defaultValue;
   std::cerr <<  std::endl << "\t\t" << "description: "<<m_description;
 };
 
 virtual void Analyze(CmdLine * cmdLine_)
 {
   m_activate = cmdLine_->Get(m_token,m_value);
 };
 
};



class Options
{
protected:
  std::list<OptionBase*> m_options;
public:
  Options(){};
  ~Options(){
    for (typename std::list<OptionBase*>::iterator it = m_options.begin(); it != m_options.end(); it++)
      {
	OptionBase * option = *it;
	delete option;
      }    
  };
  void Add(OptionBase*option_)
  {
    m_options.push_front(option_);
  };
  
  void Analyze(CmdLine * cmdLine_)
  {
    for (typename std::list<OptionBase*>::iterator it = m_options.begin(); it != m_options.end(); it++)
      {
	OptionBase * option = *it;
	option->Analyze(cmdLine_);
      }
  };

  void CheckWarning()const
  {
    for (typename std::list<OptionBase*>::const_iterator it = m_options.begin(); it != m_options.end(); it++)
      {
	const OptionBase * option = *it;
	option->CheckWarning();
      }
  };

  void CheckError()const
  {
    for (typename std::list<OptionBase*>::const_iterator it = m_options.begin(); it != m_options.end(); it++)
      {
	const OptionBase * option = *it;
	option->CheckError();
      }
  };

  void Usage() const
  {
    for (typename std::list<OptionBase*>::const_iterator it = m_options.begin(); it != m_options.end(); it++)
      {
	const OptionBase * option = *it;
	std::cerr << "\t";
	option->Usage();
	std::cerr << std::endl;
      }
  };

};


class Program
{
protected:
  Options	m_options;
  Options	m_commonOptions;
  CmdLine*	m_cmd;
  pid_t		m_pid; 		/*!< program id */
  int		m_tid; 		/*!< thread  id */
  std::string	m_progname;	/*!< name of the program */
  bool		m_verbose;
  bool		m_help;
  std::string 	m_ofilename;
public:
  
  virtual void Usage() const 
  {    
    std::cerr << m_progname << " " << std::endl;
    m_options.Usage();
    m_commonOptions.Usage();
  };
  virtual void Main() = 0;

  virtual bool HasVerbose() const
  {
    return this->m_verbose;
  };

  int GetNumInputFiles() const { return this->m_cmd->GetNumArgs()-1; };
  std::string GetInputFilename(int inputFileIndex_) const { return this->m_cmd->GetArg(inputFileIndex_+1); };

  void Info() const
  {
    if (this->HasVerbose())
      {
	Logging::Message::Get() << "Command:" << *m_cmd << std::endl;
	Logging::Message::Get() << "Plateform:" << Mpi::GetProcessorName() << std::endl;
	Logging::Message::Get() << "Version:" << Mpi::GetVersion() << std::endl;
	Logging::Message::Get() << "NumCpus:" << Mpi::GetNumCpus() << std::endl;
      }
  };

  std::string GetOfilename()const{return m_ofilename;};
  void AddOption(OptionBase*option_)
  {
    m_options.Add(option_);
  };

 Program(int argc_,char * argv_[],bool outputOptionRequired_=false) :
    m_pid(getpid()),
    m_tid(0),
    m_progname(argv_[0])
      {	
	Mpi::Init(&argc_,&argv_);
	Logging::Message::Create(argv_[0]);
	Logging::Warning::Create(argv_[0]);
	Logging::Error::Create(argv_[0]);
	if (Mpi::IsMaster())
	  {
	    m_cmd = new CmdLine(argc_,argv_);
	    m_commonOptions.Add(new Option<bool>(&m_help,"-h",false,true,"Get help."));
	    m_commonOptions.Add(new Option<bool>(&m_verbose,"-v",false,true,"Get verbose."));
	    m_commonOptions.Add(new Option<std::string>(&m_ofilename,"-o",std::string(argv_[0])+std::string(".out"),!outputOptionRequired_,"Get output filename."));
	  }
	Mpi::Barrier();
	// INITIALIZE MPI
      };
  
  virtual ~Program()
    {

      Mpi::Barrier();
      if (Mpi::IsMaster())
	{
	  if (this->m_cmd!=NULL)
	    {
	      delete this->m_cmd;
	    }
	}
      Mpi::Finalize();
      Multiton<Logging::Base,int>::KillInstances();
    };

  void Run()
  {    
    if (Mpi::IsMaster())
      {
	m_commonOptions.Analyze(m_cmd);      
	if (m_help)
	  {
	    this->Usage();
	    exit(1);
	  }
	m_options.Analyze(m_cmd);      
	
	m_options.CheckError();
	m_commonOptions.CheckError();
	m_options.CheckWarning();
	m_commonOptions.CheckWarning();
	this->m_cmd->CheckInvalid();
	if (this->HasVerbose())
	  {
	    this->Info();
	  }
      }
    Mpi::Barrier();
    this->Main();
  };

 void Msg(const char* msg_,...)
 {
   if (true==this->m_verbose)
     {
       char msg[256];
       { va_list args;
	 va_start(args,msg_);
	 vsprintf(msg,msg_,args);
	 va_end(args); }
       std::cout << "// "<< this->m_progname <<"." << this->m_pid << ": " << msg   << std::endl;
     }
 };

};

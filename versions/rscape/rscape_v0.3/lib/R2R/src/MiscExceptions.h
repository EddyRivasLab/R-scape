/*
This file copyright (c) 2009-2012, Zasha Weinberg
All rights reserved.

This copyrighted source code is freely 
distributed under the terms of the GNU
General Public License.  See the file
LICENSE in this directory for details.
*/
/*
Miscellaneous exception classes derived from std::exception

  Requires:
  <exception>,<string>  (STD C++)
*/

#ifndef MISC_EXCEPTIONS_INCLUDED
#define MISC_EXCEPTIONS_INCLUDED

// this function should be somewhere else, but it's here for historical reasons
// this function is not so efficient
std::string stringprintf(const char *format,...);

// stores the message as a string
class SimpleStringException : public std::exception {
protected:
	std::string msg;

	// NOTE: you're not supposed to throw pointers to these exceptions (unlike my previous
	// exception classes), so use of operator new probably means that I'm being absent-minded
	// thus, I'm making this non-public
protected: void * operator new (size_t bytes);

public:
	SimpleStringException (const char *format,...);
	SimpleStringException (const std::string& s);
	SimpleStringException (const SimpleStringException& t);
	~SimpleStringException () throw ();

    const char *what() const throw();
};

// 'fopen' call failed
class FopenException : public SimpleStringException {
protected:
	static std::string BuildErrorMessage(const char *fileName);
public:
	FopenException (const char *fileName);
	~FopenException () throw ();
};
// convenience: throws exception on failure
extern FILE *ThrowingFopen (const char *fileName,const char *mode);

// some ANSI C call failed
class ANSICLibException : public SimpleStringException {
protected:
	static std::string BuildErrorMessage(const char *description,const char *failedFunctionName);
public:
	ANSICLibException (const char *description,const char *failedFunctionName);
	~ANSICLibException () throw ();
};

// convenience function for ANSI-C errors, using the errno/strerror interface
std::string GetAnsiCErrorMessage (void);
// same, for WIN32 errors, using GetLastError
std::string GetWin32ErrorMessage (void);

// assert, even in release mode; throw exception on failure
#ifdef _DEBUG
// in debug mode, use regular assert
#define assertr(exp) assert(exp); if (!(exp)) { throw SimpleStringException("Internal error (release mode assertion failed \"%s\") %s:%d",#exp,__FILE__,__LINE__); } // throwing in _DEBUG mode means we don't have to bother with a 'return' statement -- a slight convenience
#else
#define assertr(exp) if (!(exp)) { throw SimpleStringException("Internal error (release mode assertion failed \"%s\") %s:%d",#exp,__FILE__,__LINE__); }
#endif

#define NotImplemented() 	throw SimpleStringException("Not implemented at %s:%d",__FILE__,__LINE__);
#define InternalError() 	throw SimpleStringException("Internal error at %s:%d",__FILE__,__LINE__);

#endif

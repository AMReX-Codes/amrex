//BL_COPYRIGHT_NOTICE

//
// $Id: Tracer.H,v 1.2 1997-09-17 22:51:21 lijewski Exp $
//

#ifndef BL_TRACER_H
#define BL_TRACER_H

//
//@Man:
//@Memo: Trace program execution.
/*@Doc:
        The definition of a trace class used for monitoring program and
        function execution.  It is meant to be used in a very stylized manner.
        By setting envirionment variables one can enable and disable complete
        program tracing, or the tracing of selected functions.
*/

class Tracer
{
  public:

    /*@ManDoc:
      Constructor.  The argument to the first Tracer() object in the
      execution path of your program is taken to be the name of your
      program.  If an environment variable of that program name, with
      "TRACER" appended to it, exists, then program tracing is turned
      on.  If that environment variable has an argument, that value
      is taken to be a string, with only function names containing
      that substring (ala strstr()) being traced.
    */
    Tracer (const char* function);
    //
    //@ManDoc: Destructor.
    //
    ~Tracer ();

  private:

    void init ();

    const char*        m_func;
    static bool        m_trace;
    static int         m_count;
    static const char* m_program;
    static char*       m_substr;
};

//
// Macro Interface to Tracer Object.
//
// Instead of using a Tracer object directly, we allow people to use a TRACER
// object.  This is a macro which can be completely removed from source code
// by compiling with -NDEBUG.
//

#ifndef NDEBUG
#define TRACER(s) Tracer  _trace(s)
#else
#define TRACER(s)
#endif

#endif /*BL_TRACER_H*/

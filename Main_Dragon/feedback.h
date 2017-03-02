//----------------------------------------
//header necessary to handle feedback information
//By Chunhua Liao, 3/5/2004
//
// Extracted from several headers in open64
//-------------------------------------------

#ifndef fb_INCLUDED
#define fb_INCLUDED


typedef signed int	INT;	/* The natural integer on the host */
typedef signed int	INT8;	/* Use the natural integer */
typedef signed int	INT16;	/* Use the natural integer */
typedef signed int	INT32;	/* The natural integer matches */
typedef signed long long INT64;	
typedef unsigned int	UINT;	/* The natural integer on the host */
typedef unsigned int	UINT8;	/* Use the natural integer */
typedef unsigned int	UINT16;	/* Use the natural integer */
typedef unsigned int	UINT32;	/* The natural integer matches */
typedef unsigned long long UINT64;
typedef int		BOOL;	/* Natural size Boolean value */
typedef signed char	mINT8;	/* Avoid - often very inefficient */
typedef signed short	mINT16;	/* Use a 16-bit integer */
typedef signed int	mINT32;	/* The natural integer matches */
typedef signed long long mINT64;
typedef unsigned char	mUINT8;	/* Use the natural integer */
typedef unsigned short	mUINT16;/* Use a 16-bit integer */
typedef unsigned int	mUINT32;/* The natural integer matches */
typedef unsigned long long mUINT64;
typedef unsigned char	mBOOL;	/* Minimal size Boolean value */

/* Define the IDTYPE used by wopt */
typedef mUINT32 IDTYPE;


//-------------------------------------------------------------------

enum FB_EDGE_TYPE {
  FB_EDGE_UNINIT           =  0,  // NOTE: If any of these are changed, then
  FB_EDGE_INCOMING         =  1,  // the array FB_EDGE_NAMES in fb_info.cxx
  FB_EDGE_OUTGOING         =  2,  // must also be updated!
  FB_EDGE_ENTRY_OUTGOING   =  3,
  FB_EDGE_BRANCH_TAKEN     =  4,
  FB_EDGE_BRANCH_NOT_TAKEN =  5,
  FB_EDGE_LOOP_ZERO        =  6,
  FB_EDGE_LOOP_POSITIVE    =  7,
  FB_EDGE_LOOP_OUT         =  8,
  FB_EDGE_LOOP_BACK        =  9,
  FB_EDGE_LOOP_EXIT        = 10,  // EXIT    == ZERO     + OUT
  FB_EDGE_LOOP_ITERATE     = 11,  // ITERATE == POSITIVE + BACK
  FB_EDGE_CIRCUIT_LEFT     = 12,
  FB_EDGE_CIRCUIT_RIGHT    = 13,
  FB_EDGE_CIRCUIT_NEITHER  = 14,
  FB_EDGE_CALL_INCOMING    = 15,
  FB_EDGE_CALL_OUTGOING    = 16,
  FB_EDGE_CALL_INOUTSAME   = 17,
  FB_EDGE_IO_OUTGOING      = 18,
  FB_EDGE_IO_ESCAPE_BASE   = 19,
  FB_EDGE_SWITCH_DEFAULT   = 22,  // 19 + FB_IO_ESCAPE_EDGES_MAX
  FB_EDGE_SWITCH_BASE      = 23   // must be last
};
extern const char *FB_EDGE_NAMES[];

// ====================================================================
// FB_FREQ_TYPE - frequency types
// ====================================================================

// The FB_FREQ_TYPE integer values are chosen to satisfy two rules:
// (1) When two FB_FREQs are combined, the resulting FB_FREQ has normally
//     FB_FREQ_TYPE equal to the minimum of those of the two FB_FREQ
//     addends; and
// (2) EXACTness and KNOWNness can be ascertained by comparing the
//     FB_FREQ_TYPE to zero.

enum FB_FREQ_TYPE {
  FB_FREQ_TYPE_EXACT   =  1,
  FB_FREQ_TYPE_GUESS   =  0,
  FB_FREQ_TYPE_UNKNOWN = -1,
  FB_FREQ_TYPE_UNINIT  = -2,
  FB_FREQ_TYPE_ERROR   = -3
};

#define FGEDGE_TYPE_NAME_LENGTH 20  // buffer length required for
				     //   FB_EDGE_TYPE_sprintf
#define FB_EDGE_TYPE_NAME_LENGTH 20

extern INT  FB_EDGE_TYPE_sprintf( char *buffer, const FB_EDGE_TYPE fb_type );



class FB_FREQ {
public:

  FB_FREQ_TYPE  _type;
  float         _value;

  
  FB_FREQ( FB_FREQ_TYPE type, float value )
    : _type( type ), _value( value ) { }

public:

  // Constructor methods

  FB_FREQ()
    : _type( FB_FREQ_TYPE_UNINIT ),
      _value( (float) FB_FREQ_TYPE_UNINIT ) {}
  /*
  FB_FREQ( float value, bool exact )
    : _type( exact ? FB_FREQ_TYPE_EXACT : FB_FREQ_TYPE_GUESS ),
      _value( value )
    { Is_True( value >= 0.0, ( "FB_FREQ: negative value %f", value ) ); }

  FB_FREQ( INT64 value )
    : _type( FB_FREQ_TYPE_EXACT ),
      _value( (float) value )
    { Is_True( value >= 0, ( "FB_FREQ: negative value %lld", value ) ); }
  */
  FB_FREQ( FB_FREQ_TYPE type )
    : _type( type ),
      _value( type >= 0 ? 0.0 : (float) type ) {
    // Is_True( FB_FREQ_TYPE_IS_VALID( type ),
    //   ( "FB_FREQ: invalid type %d", type ) );
  }

 // Printing methods

  void Print( FILE *fp ) const {
    switch ( _type ) {
    case FB_FREQ_TYPE_EXACT:
      fprintf( fp, "%g!", _value );
      break;
    case FB_FREQ_TYPE_GUESS:
      fprintf( fp, "%g?", _value );
      break;
    case FB_FREQ_TYPE_UNKNOWN:
      fprintf( fp, "unknown" );
      break;
    case FB_FREQ_TYPE_UNINIT:
      fprintf( fp, "uninitialized" );
      break;
    case FB_FREQ_TYPE_ERROR:
      fprintf( fp, "error" );
      break;
    default:
      // Is_True( FALSE, ("FB_FREQ: Unexpected type %d", _type ));
      printf("Fatal!: FB_FREQ: Unexpected type %d", _type);
      break;
    }
  }

  void Print_simple (FILE * fp) const {
  	fprintf(fp, "_type = %d   |  _value = %f \n", _type, _value);
  	}

  INT Sprintf( char *buffer ) const {
    INT length = 0;
    switch ( _type ) {
    case FB_FREQ_TYPE_EXACT:
      length = sprintf( buffer, "%g!", _value );
      break;
    case FB_FREQ_TYPE_GUESS:
      length = sprintf( buffer, "%g?", _value );
      break;
    case FB_FREQ_TYPE_UNKNOWN:
      length = sprintf( buffer, "unknown" );
      break;
    case FB_FREQ_TYPE_UNINIT:
      length = sprintf( buffer, "uninitialized" );
      break;
    case FB_FREQ_TYPE_ERROR:
      length = sprintf( buffer, "error" );
      break;
    default:
      // Is_True( FALSE, ("FB_FREQ: Unexpected type %d", _type ));
      printf("FB_FREQ: Unexpected type %d", _type);
      break;
    }
    return length;
  }


};


// Some FB_FREQ constants.  For unknown and uninitialized frequencies,
// use these instead of invoking an FB_FREQ constructor.

//const FB_FREQ FB_FREQ_ZERO(    0.0, true /*EXACT*/  );
//const FB_FREQ FB_FREQ_UNKNOWN( FB_FREQ_TYPE_UNKNOWN );
const FB_FREQ FB_FREQ_UNINIT(  FB_FREQ_TYPE_UNINIT  );
//const FB_FREQ FB_FREQ_ERROR(   FB_FREQ_TYPE_ERROR   );


#endif  /*  fb_INCLUDED */

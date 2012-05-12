// Include file for the R3 Track

// Track defintion

typedef enum {
	TRACK_STRAIGHT,
	TRACK_APPROACH_RIGHT,
	TRACK_APPROACH_LEFT,
	TRACK_EXIT_RIGHT,
	TRACK_EXIT_LEFT,
	TRACK_TURN_RIGHT,
	TRACK_TURN_LEFT,
	NUM_TRACK_TYPES
} R3TrackType;

struct R3Track {
  public:
    // Constructor functions
    R3Track(void);
  
    // I/O functions
    int Read(const char *filename);
  
  public:
	bool isCovered;			// tunnel?
	R3TrackType type;		// straight, turn, or approach
    R3Point start;			// position of beginning of along vector
	R3Point end;			// position of end of along vector
    R3Point center;         // center of curved path
	R3Vector startNormal;	// normal at beginning
	R3Vector endNormal;		// normal at end
    R3Vector along;			// in world coordinates
    R3Plane endPlane;		// end of track segment
	R3Vector side;			// vector from center to right edge
	double cof;				// coefficient of friction
    double R;               // radius full curved path
    double r;               // radius of the path itself
};

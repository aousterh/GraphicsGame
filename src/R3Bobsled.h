
// Include file for the R3 Bobsled


// Bobsled defintion

struct R3Bobsled {
  public:
    // Constructor functions
    R3Bobsled(void);
  
    // I/O functions
    int Read(const char *filename);
  
  public:
    R3Point position;   // position of center of volume of bbox
    R3Vector velocity;  // in world coordinates
    double mass;
    R3Material *material;
    R3Mesh *mesh;
    R3Track *track;     // current section of track that this bobsled is on
    R3Camera *camera;     // bobsled's camera
  
    // NOTE: transform camera position whenever you transform bobsled position!
    // scene graph, bobsled position, camera position
};

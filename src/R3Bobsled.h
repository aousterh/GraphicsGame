// Include file for the R3 Bobsled

// Bobsled defintion
struct R3Bobsled {
public:
  // Constructor functions
  R3Bobsled(void);
  
  // Update bobsled position
  void UpdateBobsled(R3Node *node, double current_time, double delta_time,
                     bool force_left, bool force_right);
  
  // I/O functions
  int Read(const char *filename);
  
public:
  R3Point position;   // position of center of volume of bbox
  R3Vector velocity;  // in world coordinates
  double mass;
  R3Node *node;        // node that stores the meshes
  //  R3Track *track;     // current section of track that this bobsled is on
  
  // update scene graph, bobsled position when bobsled moves
};
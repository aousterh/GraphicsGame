
<<<<<<< HEAD
// Include file for the R3 Track


// Track defintion

struct R3Track {
  public:
    // Constructor functions
    R3Track(void);
  
        // I/O functions
    int Read(const char *filename);
  
  public:
    R3Point position;   // position of center of volume of bbox
    R3Vector velocity;  // in world coordinates
    double mass;
    R3Material *material;
    R3Mesh *helmets;
    R3Mesh *masks;
    
    R3Track *track;     // current section of track that this bobsled is on
    R3Camera *camera;     // bobsled's camera
    
    double main_theta;
    double inside_theta;
  
    // NOTE: transform camera position whenever you transform bobsled position!
    // scene graph, bobsled position, camera position
};
=======
// Update bobsled position
void UpdateBobsled(R3Scene *scene, double current_time, double delta_time, bool force_left, bool force_right);
>>>>>>> 1326cd87873ec13084d934af14a7f58bf700b7f2

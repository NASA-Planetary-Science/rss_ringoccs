
/*  Data types for two and three dimensional points, respectively.            */
typedef struct rssringoccs_TwoVector {
    double dat[2];
} two_vec;

typedef struct rssringoccs_ThreeVector {
    double dat[3];
} rssringoccs_ThreeVector;

extern rssringoccs_TwoVector
rssringoccs_Normalize_TwoVector(rssringoccs_TwoVector P);

extern rssringoccs_ThreeVector
rssringoccs_Normalize_ThreeVector(rssringoccs_ThreeVector P);

extern rssringoccs_ThreeVector
rssringoccs_Cross_Product(rssringoccs_ThreeVector P, rssringoccs_ThreeVector Q);

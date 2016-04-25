//==============================================================================
#ifndef BEAM_HPP
#define BEAM_HPP

#ifndef __LIST_H
#include<list>
#endif

#ifndef GEOMETRY_HPP
#include "Geometry.hpp"
#endif

//extern double d;


//==============================================================================
// 3D polygon + Johns matrix
/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** 
@brief The Beam class consist of 3D polygon, propagation direction and it's history, Jones matrix and so on.
*/
class Beam {
private:
	matrixC mt;             ///< Jones matrix of the beam
public:
	std::list<Point3D> v;   ///< The beam's vertexes
	std::list<unsigned int> path; ///<beam's trajectory

	double	lng,           ///< optical path of the beam
			D;             ///< current position of phase front from Ax+By+Cz+D=0
	Point3D	e,             ///< the basis of polarization plane
			r,             ///< direction of the beam in 3D space
			T,F,N; 

	bool isInternal = true; ///< internal or external

	Beam(void) : mt(2,2), lng(0) 
		{ this->mt.Identity(); }
	Beam(const std::list<Point3D>& _v, double l = 0) :
		mt(2,2), v(_v), lng(l)
		{ this->mt.Identity(); }
	~Beam(void) {}
	Beam&  operator=(const Beam& b)
	{ 
		this->mt = b.mt;
		this->v = std::move(b.v);
		this->path = std::move(b.path);

		this->lng = b.lng;
		this->D = std::move(b.D);
		this->e = std::move(b.e);
		this->r = std::move(b.r);
		this->T = std::move(b.T);
		this->F = std::move(b.F);
		this->N = std::move(b.N);
		this->isInternal = b.isInternal;
		return *this;
	}

	complex   operator() (unsigned int n, unsigned int m) const ///< Returns Jones matrix element of the beam
		{ return this->mt[n][m]; }
	complex&  operator() (unsigned int n, unsigned int m)	///< Return reference to Jones matrix element
		{ return this->mt[n][m]; }
	matrixC   operator() (void) const  { return this->mt; } ///< Returns Jones matrix of the beam
	matrixC&  operator() (void)        { return this->mt; } ///< Return reference to Jones matrix 
	void      Clear(void)              { this->v.clear(); } ///< Clears beam vertex list
	std::list<Point3D>::const_iterator  Begin(void) const 
		{ return this->v.begin(); }
	std::list<Point3D>::const_iterator  End(void) const
		{ return this->v.end(); }
	std::list<unsigned int>::const_iterator  BeginP(void) const
		{ return this->path.begin(); }
	std::list<unsigned int>::const_iterator  EndP(void) const
		{ return this->path.end(); }

	Polyg     Projection(const double * const, int) const;
	Beam&     Rotate(const Point3D &k, const Point3D &Ey);  ///< Rotate Jones matrix of the beam to basis of scattering plane
	Beam      RotatePlane(const Point3D & NewE); ///< Rotate Jones matrix to new basis
	void      PushFront(const Point3D& p) { this->v.push_front(p); }
	void      PushBack(const Point3D& p)  { this->v.push_back(p); }
	void      PushFrontP(const int& f) { this->path.push_front(f); }
	void      PushBackP(const int& f)  { this->path.push_back(f); }	
	SphCrd    Spherical(void) const;
	// friends
	friend unsigned int  Size(const Beam& b)  { return b.v.size(); }
	friend unsigned int  SizeP(const Beam& b) { return b.path.size(); }
	friend double&  OpticalPath(Beam& b) { return b.lng; } ///<returns optical path
	friend Point3D  CenterOfBeam(const Beam &); ///< returns center of the beam
	friend double   CrossSection(const Beam &); 
	friend double   AreaOfBeam(const Beam &);
	friend void     OutBeam(char*, const Beam &);
	friend Point2D  Proj_Vertex(const Beam &, Point3D &);
	friend void Rot(const Point3D&, const Point3D&, const Point3D&, matrixC &); ///< Rotates coordinates system
};

/** @}
*/
#endif












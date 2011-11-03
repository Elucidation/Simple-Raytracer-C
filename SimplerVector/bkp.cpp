#include <stdio.h>
//#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <ctime>
using namespace std;

#define ResolutionX 200 // Damnit.
#define ResolutionY 200 // Keep equal values, not sure how to do diff sized X & Y yet -_-
#define LIGHTS 1
#define SPHERES 0
#define SLABS 64*64

/////////////////////////////////////////////////////////////////////////////////////////
template <class T>
class Stack {
public:
	Stack(int sz);
	~Stack();
	void Push(T value);
	bool Full();
	void show();
//private:
	int size;
	int top;
	T *stack;
};

template <class T>
Stack<T>::Stack(int sz)
{
	size = sz;
	top = 0;
	stack = new T[size];
}

template <class T>
Stack<T>::~Stack()
{
	delete [] stack;
}

template <class T>
void Stack<T>::Push(T value)
{
	assert(!Full());
	stack[top++] = value;
}

template <class T>
bool Stack<T>::Full()
{
	return (top == size);
}

/////////////////////////////////////////////////////////////////////////////////////////



struct V
{
	float x;
	float y;
	float z;
};

struct Color
{
	float r;
	float g;
	float b;
};

struct Prim
{
	char type;
	int index;
};

struct Sphere
{
	V center;
	float radius;
	Color color;
	float transparant;
	float reflect;
};

struct Slab
{
	V minB;
	V maxB;
	Color color;
	float transparant;
	float reflect;
};

struct Ray
{
	V origin;
	V dir;
};

struct Cam
{
	V origin;
	V dir;
	V* rays;//= new V[ResolutionX*ResolutionY];
};

struct Light
{
	V pos;
	Color col;
};
struct ColorMap
{
	Color* colors;// = new Color[ResolutionX*ResolutionY];
};


void setC(Color &color,float r, float g, float b)
{
	color.r = r;
	color.g = g;
	color.b = b;
}
void setV(V &vector, float const x, float const y, float const z)
{
	vector.x = x;
	vector.y = y;
	vector.z = z;
}
void showV(V const vector)
{
	cout << "<" << vector.x << "," << vector.y << "," << vector.z << ">" << endl;
}
void addV(V &out, V const a, V const b)
{
	out.x = a.x + b.x;
	out.y = a.y + b.y;
	out.z = a.z + b.z;
}
void minV(V &out, V const a, V const b)
{
	out.x = a.x - b.x;
	out.y = a.y - b.y;
	out.z = a.z - b.z;
}
void dotV(float &out, V const a, V const b)
{
	out = a.x*b.x + a.y*b.y + a.z*b.z;
}
void crossV(V &out, V const a, V const b)
{
	// a b c // x y z
	// d e f // x y z
	// i = b*f - c*e
	// j = -(a*f - c*d)
	// k = a*e - b*d
	out.x = a.y*b.z - a.z*b.y;
	out.y = a.z*b.x - a.x*b.z;
	out.z = a.x*b.y - a.y*b.x;
}
void magV(float &out, V const a)
{
	out = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

void normV(V &out, V const a)
{
	float mag;
	magV(mag,a);
	out.x = a.x / mag;
	out.y = a.y / mag;
	out.z = a.z / mag;
}

void multStoV(V &out, V const a, float const b)
{
	out.x = a.x * b;
	out.y = a.y * b;
	out.z = a.z * b;
}

void compBinA(float &out, V const a, V const b)
{
	// comp = a*b / |a|
	dotV(out,a,b); // a*b
	float mag;
	magV(mag,a); // |a|
	out /= mag; // a*b / |a|
}
void projBonA(V &out, V const a, V const b)
{
	float comp;
	compBinA(comp,a,b);
	V dir;
	normV(dir,a);
	multStoV(out,dir,comp);
}

float intersectRaySlab(V const rayOrigin, V const rayDir, V const slabMin, V const slabMax)
{
	float Tnear = -1000000;
	float Tfar = 1000000;

	float _l[3];
	float _h[3];
	float _o[3];
	float _d[3];

	_l[0] = slabMin.x;
	_l[1] = slabMin.y;
	_l[2] = slabMin.z;

	_h[0] = slabMax.x;
	_h[1] = slabMax.y;
	_h[2] = slabMax.z;

	_o[0] = rayOrigin.x;
	_o[1] = rayOrigin.y;
	_o[2] = rayOrigin.z;

	_d[0] = rayDir.x;
	_d[1] = rayDir.y;
	_d[2] = rayDir.z;

	for(int i=0; i<3; i++)
	{
		if (_d[i] == 0) // Parallel to X plane?
		{
			if (_o[i] < _l[i] || _o[i] > _h[i])
			{return -1;} // False if not in the area
		}
		else // Not Parallel
		{
			float t1 = (_l[i] - _o[i]) / _d[i];
			float t2 = (_h[i] - _o[i]) / _d[i];
			if (t1 > t2)
			{
				float ot1 = t1; // holder float
				t1 = t2;
				t2 = ot1;
			}
			if (t1 > Tnear) { Tnear = t1; }
			if (t2 < Tfar) { Tfar = t2; }
			if (Tnear > Tfar) { return -1; } // box is missed
			if (Tfar < 0) { return -1; } // box behind ray miss
		}
	}
	return Tnear; // Tfar is exit point, Tnear is entrance
}

float simpleIntersectRaySlab(Ray const r, Slab const s)
{
	return intersectRaySlab(r.origin,r.dir,s.minB,s.maxB);
}

void getSlabNormal(V &slabNorm,V const P,Slab const obj)
{
	if (abs(obj.minB.x - P.x) < 0.00001)
	{
		setV(slabNorm,-1,0,0); // Front (Facing towards camera if camera facing in x axis direction)
	}
	else if (abs(obj.maxB.x - P.x) < 0.00001)
	{
		setV(slabNorm,1,0,0); // Back
	}
	else if (abs(obj.minB.y - P.y) < 0.00001)
	{
		setV(slabNorm,0,-1,0); // Left
	}
	else if (abs(obj.maxB.y - P.y) < 0.00001)
	{
		setV(slabNorm,0,1,0); // Right
	}
	else if (abs(obj.minB.z - P.z) < 0.00001)
	{
		setV(slabNorm,0,0,-1); // Bottom
	}
	else if (abs(obj.maxB.z - P.z) < 0.00001)
	{
		setV(slabNorm,0,0,1); // Top
	}
	else
	{
		setV(slabNorm,1,0,0); // Bad Position, Choose something innocuous like Back
	}
}
float intersectRaySphere(V const rayOrigin, V const rayDir, V const sphereCenter, float const sphereRadius)
{
	V dst;
	minV(dst,rayOrigin,sphereCenter);
	float B;
	float C;
	float D;
	dotV(B,dst,rayDir);
	dotV(C,dst,dst);
	C -= sphereRadius*sphereRadius;
	D = B*B - C;
	return D > 0 ? -B - sqrt(D) : -1; // return t of intersect = r0 + t*rd; else return -1 for false no intersection
}

float simpleIntersectRaySphere(Ray const r, Sphere const s)
{
	return intersectRaySphere(r.origin,r.dir,s.center,s.radius);
}
void getSphereNormal(V &norm,V const intersect ,Sphere const obj)
{
	minV(norm,intersect,obj.center);
	normV(norm,norm);
}

void makeCameraVecs(Cam &c,float focalLength)
{
	c.rays = new V[ResolutionX*ResolutionY];
	//float focalLength = 10;
	V Cview = c.dir;
	//minV(Cview,c.dir,c.origin); // This was if c.dir was endpt so endpt-origin gets ray (and norm it)
	V initTopDir;
	setV(initTopDir,0,0,1);
	V HorizontalDir;
	crossV(HorizontalDir,Cview,initTopDir);
	V TopDir;
	crossV(TopDir,Cview,HorizontalDir);

	// Forward dir = Cview
	// Horizontal dir = HorizontalDir
	// Vertical dir = TopDir;

	V Hstep;
	V Vstep;
	normV(Hstep,HorizontalDir); // Not sure if crossing 2 normals gets a normal, just in case.
	normV(Vstep,TopDir);  // Ditto
	multStoV(Hstep,Hstep,10.0/ResolutionX);
	multStoV(Vstep,Vstep,10.0/ResolutionY);

	V endPt;
	multStoV(endPt,c.dir,focalLength);
	addV(endPt,c.origin,endPt); // Setting up point of middle of the screen

	for (int i=0;i<ResolutionX;i++)
	{
		for (int j=0;j<ResolutionY;j++)
		{
			V cRay;
			V Hdiff;
			V Vdiff;
			multStoV(Hdiff,Hstep,i-ResolutionX/2); // Hdiff = Hstep * (i-ResolutionX/2)
			multStoV(Vdiff,Vstep,j-ResolutionY/2); // Vdiff = Vstep * (j-ResolutionY/2) // Half to make endpt the mid
			addV(cRay,endPt,Hdiff); // cRay = endPt + Hdiff
			addV(cRay,cRay,Vdiff); // cRay += Vdiff
			//V** a = c.rays[i][j]; // Can't figure out a way to access vectors
			//V* b = *a;
			//b->x; // This works
			//cout << b.x << endl; // This crashes
			//V a = **c.rays[i][j]; // This does work, wtfffff
			//showV(a); // This doesn't work
			//c.rays[0].x = 5; // Huzzah it works now //Wtf this crashes too -_-

			c.rays[i*ResolutionX + j] = cRay; // Make THIS work somehow
		}
	}

}

void getReflectionRay(V &reflectRay,V const ndir,V const objNormAtPt)
{
	V cVec;
	V nVec;
	projBonA(nVec,objNormAtPt,ndir); // Get Normal Force
	minV(cVec,ndir,nVec); // Get Crossing Force
	multStoV(reflectRay,nVec,-1.0);
	addV(reflectRay,reflectRay,cVec);
}



Stack<Light> LightList(LIGHTS);
Stack<Sphere> SphereList(SPHERES);
Stack<Slab> SlabList(SLABS);

float shadowRaytrace(V const point,Light const clight)
{
	V sDir;
	float lShadow = 1; // Shadow from this light
	minV(sDir,clight.pos,point); // Setup Ray from point to light
	float dist;
	magV(dist,sDir);
	normV(sDir,sDir);
	// Slabs
	for (int i=0;i<SlabList.top;i++)
	{
		float intersect = intersectRaySlab(point, sDir, SlabList.stack[i].minB, SlabList.stack[i].maxB);
		if (intersect > 0.0000001 && intersect < dist)
		{
			lShadow = 0.2;
			goto skip;
		}
	}
	// Spheres
	for (int i=0;i<SphereList.top;i++)
	{
		float intersect = intersectRaySphere(point, sDir, SphereList.stack[i].center, SphereList.stack[i].radius);
		if (intersect > 0.0000001 && intersect < dist)
		{
			lShadow = 0.2;
			goto skip;
		}
	}
	skip:;
	return lShadow;
}

///
void raytrace(Color &out, V const origin, V const ndir,Prim skip)
{
	float t = 1000000;
	Color c;
	setC(c,0.2,0.2,0.2);
	V intersectionPt; // P = P0 + t * Pd
	V objNormAtPt;
	Color objColor;
	float transparant = 0;
	float reflect = 0;
	V exitPt;
	Prim obj;
	obj.index = -1;
	obj.type = 'n';
	for (int i=0;i<SphereList.top;i++)
	{
		float intersect = intersectRaySphere(origin, ndir, SphereList.stack[i].center, SphereList.stack[i].radius);
		if (intersect >= 0 && intersect < t)
		{
			if (skip.type == 's' && skip.index == i) {}
			else{
				t = intersect;
				multStoV(intersectionPt,ndir,t); // P = Pd * t
				addV(intersectionPt,intersectionPt,origin); // P = P + P0
				getSphereNormal(objNormAtPt,intersectionPt,SphereList.stack[i]);
				V bumpV;
				//minV(bumpV,intersectionPt,SphereList.stack[i].center); // bump = P - SphereCenter
				//normV(objNormAtPt,objNormAtPt);
				multStoV(bumpV,objNormAtPt,0.000001); // bump = norm * bumpFloat
				addV(intersectionPt,intersectionPt,bumpV); // P += bump
				objColor = SphereList.stack[i].color;
				transparant = SphereList.stack[i].transparant;
				reflect = SphereList.stack[i].reflect;
				if (transparant > 0)
				{
					obj.type = 's';
					obj.index = i;
				}
			}
		}
	}
	for (int i=0;i<SlabList.top;i++)
	{
		float intersect = intersectRaySlab(origin, ndir, SlabList.stack[i].minB, SlabList.stack[i].maxB);
		if (intersect >= 0 && intersect < t)
		{
			if (skip.type == 'b' && skip.index == i) {}
			else{
				t = intersect;
				multStoV(intersectionPt,ndir,t); // P = Pd * t
				addV(intersectionPt,intersectionPt,origin); // P = P + P0
				getSlabNormal(objNormAtPt,intersectionPt,SlabList.stack[i]);
				//normV(objNormAtPt,objNormAtPt); // |bump|
				V bumpV;
				multStoV(bumpV,objNormAtPt,0.000001); // bump = bump * bumpFloat
				addV(intersectionPt,intersectionPt,bumpV); // P += bump
				//Slab obj;
				objColor = SlabList.stack[i].color;
				transparant = SlabList.stack[i].transparant;
				reflect = SlabList.stack[i].reflect;
				if (transparant > 0)
				{
					obj.type = 'b';
					obj.index = i;
				}
			}
		}
	}
	if (t > 0 && t < 1000000) // Hit Something
	{
		//V intersectionPt; // P = P0 + t * Pd
		//multStoV(intersectionPt,ndir,t - 0.001); // P = Pd * t
		//addV(intersectionPt,intersectionPt,origin); // P = P + P0
		// Check Shadow
		//out = objColor;
		setC(out,objColor.r * 0.2, objColor.g * 0.2, objColor.b * 0.2);
		for (int i=0; i< LightList.top; i++)
		{
			float diffuse;
			V lightDir;
			minV(lightDir,LightList.stack[i].pos,intersectionPt);
			normV(lightDir,lightDir);
			dotV(diffuse,lightDir,objNormAtPt);
			float shade = shadowRaytrace(intersectionPt,LightList.stack[i]);
			//Color combine;
			Color through;
			if (transparant > 0) {
				raytrace(through,intersectionPt,ndir,obj);
				//combine.r = diffuse * (1 - transparant) + through.r * transparant;
				//combine.g = diffuse * (1 - transparant) + through.g * transparant;
				//combine.b = diffuse * (1 - transparant) + through.b * transparant;
			}
			else
			{
				//setC(combine,diffuse,diffuse,diffuse);
				setC(through,0,0,0);
			}
			Color reflection;
			if (reflect > 0)
			{
				V reflectRay;
				getReflectionRay(reflectRay,ndir,objNormAtPt);
				raytrace(reflection,intersectionPt,reflectRay,obj); // Possibly make this object the skip
			}
			else
			{
				setC(reflection,0,0,0);
			}

			out.r +=  (transparant*through.r + reflect*reflection.r + LightList.stack[i].col.r * ((1-reflect)*(1-transparant)*diffuse * shade)) / LightList.top;
			out.g +=  (transparant*through.g + reflect*reflection.g + LightList.stack[i].col.g * ((1-reflect)*(1-transparant)*diffuse * shade)) / LightList.top;
			out.b +=  (transparant*through.b + reflect*reflection.b + LightList.stack[i].col.b * ((1-reflect)*(1-transparant)*diffuse * shade)) / LightList.top;
		}
		return;
	}
	else // Missed all objects
	{
		V sky;
		setV(sky,0,0,1);
		float intense;
		dotV(intense,ndir,sky);
		intense = 1 - abs(intense);
		setC(out,0.3*intense,0.3*intense,0.8*intense);
		return;
	}
}
///


void showMap(ColorMap const c)
{
	for (int i=0; i<ResolutionX; i++)
	{
		for (int j=0; j<ResolutionY; j++)
		{
			cout << c.colors[(ResolutionY - j - 1)*ResolutionX + i].r << " ";
		}
		cout << endl;
	}
}
void writeMap(ColorMap const c)
{
	ofstream rf("valuesr.txt");
	ofstream gf("valuesg.txt");
	ofstream bf("valuesb.txt");
	for (int i=0; i<ResolutionX; i++)
	{
		for (int j=0; j<ResolutionY; j++)
		{
			rf << c.colors[(ResolutionY - j - 1)*ResolutionX + i].r << " ";
			gf << c.colors[(ResolutionY - j - 1)*ResolutionX + i].g << " ";
			bf << c.colors[(ResolutionY - j - 1)*ResolutionX + i].b << " ";
		}
		rf << endl;
		gf << endl;
		bf << endl;
	}
	rf.close();
	gf.close();
	bf.close();
}



int main()
{
	// Build terrain map


	// Area 16 x 16 we shall say
	int seed = 100005;
	srand(seed);
	// Make blocks of 1x1 within it.
	cout << "Initial Seed: " << seed << ". Setting up Scene..." << endl;
	cout << LIGHTS << "Lights. " << SPHERES << " Spheres. " << SLABS << " Slabs. " << endl;
	// Lights
	Light l1;
	//          +X,+Y,+Z
	setV(l1.pos,8,8,8);
	setC(l1.col,1,1,1);
	LightList.Push(l1);

	// Spheres
	/*
	Sphere s1;
	setV(s1.center,8,0,0);
	s1.radius = 3;
	setC(s1.color,0,0,1);
	SphereList.Push(s1);
	*/

	// Slabs 
	
	Slab sl1;
	for (int i=0; i < 64; i++)
	{
		//float xpos = (float)rand()/RAND_MAX * 15.0;
		float xpos = i;
		for (int j=0; j < 64; j++)
		{
			// Ground
			//float ypos = (float)rand()/RAND_MAX * 15.0;
			float ypos = j;
			//float zpos = floor((float)rand()/RAND_MAX * 16.0);
			/*
			float zpos = 0;
			setV(sl1.minB,xpos,ypos,zpos);
			setV(sl1.maxB,xpos+1,ypos+1,zpos+1);
			//setC(sl1.color,(float)rand()/RAND_MAX*0.5+0.5,(float)rand()/RAND_MAX*0.5+0.5,(float)rand()/RAND_MAX*0.5+0.5);
			setC(sl1.color,1,1,1);
			//sl1.transparant = floor((float)rand()/RAND_MAX + 0.5) * 0.5;
			sl1.transparant = 0;
			//zpos < 3 ? sl1.reflect = 0.8 : sl1.reflect = 0;
			sl1.reflect = 1;
			SlabList.Push(sl1);
			*/

			// Floating
			float zpos = floor((float)rand()/RAND_MAX * 63.0);
			setV(sl1.minB,xpos,ypos,zpos);
			setV(sl1.maxB,xpos+1,ypos+1,zpos+1);
			setC(sl1.color,(float)rand()/RAND_MAX*0.5+0.5,(float)rand()/RAND_MAX*0.5+0.5,(float)rand()/RAND_MAX*0.5+0.5);
			sl1.transparant = 0;
			//sl1.reflect = 0;
			if (zpos < 3) {sl1.reflect = 1.0;} 
			else {sl1.reflect = 0;}
			SlabList.Push(sl1);
		}
	}

	// Camera
	Cam mainView;
	//for (int i=0;i<ResolutionX;i++){V** mainView.rays[i] = new V[ResolutionY];}
	V endPt;
	setV(endPt,32,32,10);
	setV(mainView.origin,0,-5,32);
	V cDir;
	minV(cDir,endPt,mainView.origin);
	//normV(cDir,cDir); // going to be done in 2 lines
	mainView.dir = cDir;//setV(mainView.dir,1,0,0);
	normV(mainView.dir,mainView.dir);
	cout << "Setting up Camera Rays..." << endl;
	makeCameraVecs(mainView,6);
	cout << "Done setting up Rays." << endl;
	cout << "Resolution : " << ResolutionX << " x " << ResolutionY << " pixels." << endl;

	// Intercept obj
	Prim noskip;
	noskip.index = -1;
	noskip.type = 'n';

	ColorMap pic;
	pic.colors = new Color[ResolutionX*ResolutionY];
	cout << "Beginning Raytracing scene..." << endl;
	float startTick = clock();
	float stepper = startTick;
	for (int i=0; i<ResolutionX; i++)
	{
		for (int j=0; j<ResolutionY; j++)
		{
			V rayDir;
			minV(rayDir,mainView.rays[i*ResolutionX + j],mainView.origin);
			normV(rayDir,rayDir);
			raytrace(pic.colors[i*ResolutionX + j], mainView.origin, rayDir,noskip);
		}
		if ( (clock() - stepper)/CLOCKS_PER_SEC >= 10.0 ) // If it's been 10 seconds without note
		{
			cout << "On row " << i << " / " << ResolutionX << ". " 
				 << (clock() - startTick)/CLOCKS_PER_SEC << " seconds in..." << endl;
			stepper = clock();
		}
	}
	float endTick = clock();
	delete[] mainView.rays;
	float timeToTrace = (endTick-startTick) / CLOCKS_PER_SEC;
	cout << "Done Raytracing scene. Took " << timeToTrace << " seconds." << endl;

	//showMap(pic);
	cout << "Writing information to 'values[rgb].txt'" << endl;
	writeMap(pic);
	cout << "Done writing to file 'values[rgb].txt'" << endl;
	
	delete[] pic.colors;
	
	return 0;
}
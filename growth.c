#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#pragma mark -
#pragma mark [ definitions ]

#define MAX_NEIG	32		// maximum number of neighbouring vertices that a vertex can have

typedef struct
{
	double	x,y,z;
}double3D;
typedef struct
{
	int	a,b,c;
}int3D;




typedef struct
{
	int	a,b,c,d;
}int4D;
typedef struct
{
	double	a,b,c, d,e,f, g,h,i;
}matrix;
typedef struct
{
	short	timeStamp;	// time stamp for the cell creation
	int		i;			// vertex index
	long	*next;		// next cell at the same hash
}HashCell;
typedef struct
{
	int		n;				// number of neighbours
	int		p[MAX_NEIG];	// neighbour vertex indices 
	matrix	K[MAX_NEIG];	// stiffness matrix with neighbour vertices
	matrix	A[MAX_NEIG];	// (conjugate gradient) A matrix
}Neighbours;
typedef struct
{
	int			p[4];		// vertex indices
	int			ngb[4][4];	// position in neighbours list (Each vertex has an associated
	// Neighbours structure. This matrix provides the position of
	// vertex j in the neighbours list of vertex i, for i,j=0,...,3)
	matrix		K[4][4];	// stiffness matrices
	matrix		R;			// element rotation
	double		young;		// young's modulus
	double		poisson;	// poisson's ratio
}Tetra;
typedef struct
{
	int			nt;		// number of tetrahedra
	Tetra		*t;		// tetrahedra
	
	int			np;		// number of vertices
	double3D	*p0;	// vertex coordinates at rest
	double3D	*p;		// actual vertex coordinates
	double3D	*v;		// vertex velocity
	double		*m;		// vertex mass
	double3D	*f0;	// vertex internal forces (f0)
	double3D	*fext;	// vertex external forces
	Neighbours	*ngb;	// vertex neighbours
	
	double		*fibre_length;

	HashCell	*hash;	// collision detection hash table
	int			nhash;	// number of cells allocated for collision detection
	float		h;		// cell size
	
	double3D	*b;		// (conjugate gradient) b vector
	double3D	*Ap;	// (conjugate gradient) A*p vector
	double3D	*R;		// (conjugate gradient) r
	double3D	*P;		// (conjugate gradient) p
	
	double	dt;
}Model;

#pragma mark -
#pragma mark [ linear algebra ]
double3D add3D(double3D a, double3D b)
{
	return (double3D){a.x+b.x,a.y+b.y,a.z+b.z};
}
double3D sub3D(double3D a, double3D b)
{
	return (double3D){a.x-b.x,a.y-b.y,a.z-b.z};
}
double3D sca3D(double3D a, double t)
{
	return (double3D){a.x*t,a.y*t,a.z*t};
}
double3D cro3D(double3D a, double3D b)
{
	return (double3D){a.y*b.z-b.y*a.z, b.x*a.z-a.x*b.z, a.x*b.y-b.x*a.y};
}
double dot3D(double3D a, double3D b)
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
double nor3D(double3D a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
matrix addMat(matrix a, matrix b)
{
	return (matrix) { a.a+b.a,a.b+b.b,a.c+b.c,
		a.d+b.d,a.e+b.e,a.f+b.f,
		a.g+b.g,a.h+b.h,a.i+b.i};
}
matrix trnMat(matrix a)
{
	return (matrix) { a.a, a.d, a.g,
		a.b, a.e, a.h,
		a.c, a.f, a.i};
}
double3D mulMat(matrix a, double3D b)
{
	return (double3D) {	a.a*b.x+a.b*b.y+a.c*b.z,
		a.d*b.x+a.e*b.y+a.f*b.z,
		a.g*b.x+a.h*b.y+a.i*b.z };
}
matrix scaMat(matrix a, double t)
{
	return (matrix) {	a.a*t, a.b*t, a.c*t,
		a.d*t, a.e*t, a.f*t,
		a.g*t, a.h*t, a.i*t};
}
double detMat(matrix a)
{
	return a.b*a.f*a.g + a.c*a.d*a.h + a.a*a.e*a.i - a.c*a.e*a.g - a.a*a.f*a.h - a.b*a.d*a.i;
}
matrix rkrtMat(matrix R, matrix K)
{
	return (matrix) {
		R.a*(R.a*K.a+R.b*K.d+R.c*K.g)+R.b*(R.a*K.b+R.b*K.e+R.c*K.h)+R.c*(R.a*K.c+R.b*K.f+R.c*K.i),
		R.d*(R.a*K.a+R.b*K.d+R.c*K.g)+R.e*(R.a*K.b+R.b*K.e+R.c*K.h)+R.f*(R.a*K.c+R.b*K.f+R.c*K.i),
		R.g*(R.a*K.a+R.b*K.d+R.c*K.g)+R.h*(R.a*K.b+R.b*K.e+R.c*K.h)+R.i*(R.a*K.c+R.b*K.f+R.c*K.i),
		
		R.a*(K.a*R.d+K.d*R.e+R.f*K.g)+R.b*(K.b*R.d+R.e*K.e+R.f*K.h)+R.c*(K.c*R.d+R.e*K.f+R.f*K.i),
		R.d*(K.a*R.d+K.d*R.e+R.f*K.g)+R.e*(K.b*R.d+R.e*K.e+R.f*K.h)+R.f*(K.c*R.d+R.e*K.f+R.f*K.i),
		R.g*(K.a*R.d+K.d*R.e+R.f*K.g)+R.h*(K.b*R.d+R.e*K.e+R.f*K.h)+R.i*(K.c*R.d+R.e*K.f+R.f*K.i),
		
		R.a*(K.a*R.g+K.d*R.h+K.g*R.i)+R.b*(K.b*R.g+K.e*R.h+K.h*R.i)+R.c*(K.c*R.g+K.f*R.h+R.i*K.i),
		R.d*(K.a*R.g+K.d*R.h+K.g*R.i)+R.e*(K.b*R.g+K.e*R.h+K.h*R.i)+R.f*(K.c*R.g+K.f*R.h+R.i*K.i),
		R.g*(K.a*R.g+K.d*R.h+K.g*R.i)+R.h*(K.b*R.g+K.e*R.h+K.h*R.i)+R.i*(K.c*R.g+K.f*R.h+R.i*K.i)};
}
matrix pinvxMat(matrix p, matrix x)
{
	double	det;
	double	a,b,c,d,e,f,g,h,i;
	
	det=x.b*x.f*x.g + x.c*x.d*x.h + x.a*x.e*x.i - x.c*x.e*x.g - x.a*x.f*x.h - x.b*x.d*x.i;
	
	a=(x.e*x.i - x.f*x.h);
	b=(x.c*x.h - x.b*x.i);
	c=(x.b*x.f - x.c*x.e);
	
	d=(x.f*x.g - x.d*x.i);
	e=(x.a*x.i - x.c*x.g);
	f=(x.c*x.d - x.a*x.f);
	
	g=(x.d*x.h - x.e*x.g);
	h=(x.b*x.g - x.a*x.h);
	i=(x.a*x.e - x.b*x.d);
	
	return (matrix) {
		(p.a*a + p.b*d + p.c*g)/det,
		(p.a*b + p.b*e + p.c*h)/det,
		(p.a*c + p.b*f + p.c*i)/det,
		
		(p.d*a + p.e*d + p.f*g)/det,
		(p.d*b + p.e*e + p.f*h)/det,
		(p.d*c + p.e*f + p.f*i)/det,
		
		(p.g*a + p.h*d + p.i*g)/det,
		(p.g*b + p.h*e + p.i*h)/det,
		(p.g*c + p.h*f + p.i*i)/det};
}
double3D invmxp(matrix m, double3D p)
{
	double	det;
	double	a,b,c,d,e,f,g,h,i;
	
	det=m.b*m.f*m.g + m.c*m.d*m.h + m.a*m.e*m.i - m.c*m.e*m.g - m.a*m.f*m.h - m.b*m.d*m.i;
	
	a=(m.e*m.i - m.f*m.h);
	b=(m.c*m.h - m.b*m.i);
	c=(m.b*m.f - m.c*m.e);
	
	d=(m.f*m.g - m.d*m.i);
	e=(m.a*m.i - m.c*m.g);
	f=(m.c*m.d - m.a*m.f);
	
	g=(m.d*m.h - m.e*m.g);
	h=(m.b*m.g - m.a*m.h);
	i=(m.a*m.e - m.b*m.d);
	
	return (double3D) {
		(p.x*a + p.y*d + p.z*g)/det,
		(p.x*b + p.y*e + p.z*h)/det,
		(p.x*c + p.y*f + p.z*i)/det};
}
matrix ortMat(matrix r)
{
	/*
	 double3D	x={r.a,r.d,r.g};
	 double3D	y={r.b,r.e,r.h};
	 double3D	z=cro3D(x,y);
	 double		nx=nor3D(*(double3D*)&r);
	 double		ny=nor3D(*(double3D*)&(r.d));
	 double		nz=nor3D(z);
	 
	 return (matrix) {	x.x/nx, y.x/ny, z.x/nz,
	 x.y/nx, y.y/ny, z.y/nz,
	 x.z/nx, y.z/ny, z.z/nz};
	 */
	double3D	x={r.a,r.d,r.g};
	double3D	y={r.b,r.e,r.h};
	x=sca3D(x,1/nor3D(x));
	double3D	z=cro3D(x,y);
	z=sca3D(z,1/nor3D(z));
	y=cro3D(z,x);
	y=sca3D(y,1/nor3D(y));
	
	return (matrix) {
		x.x, y.x, z.x,
		x.y, y.y, z.y,
		x.z, y.z, z.z};
}
matrix vecs2mat(double3D a, double3D b, double3D c)
{
	return (matrix) { a.x,b.x,c.x,
		a.y,b.y,c.y,
		a.z,b.z,c.z };
}
#pragma mark -
#pragma mark [ model: initialise ]
void model_setTetra(Model *m,int tetraIndex,int a,int b,int c,int d,double E,double nu,double rho)
/*
 Initialise the tetrahedron with index tetraIndex in the model m. The
 tetrahedron vertices are those with indices a,b,c,d. Young's
 coefficient E, Poisson's factor nu, and mass density rho. The tetra-
 hedron mass is computed from its volume and rho, and divided among
 the four vertices
 */
{
	int		i;
	Tetra	*t=&(m->t[tetraIndex]);
	double	V;
	
	*(int4D*)(t->p)=(int4D){a,b,c,d};
	V=detMat(vecs2mat(sub3D(m->p0[t->p[1]],m->p0[t->p[0]]),
					  sub3D(m->p0[t->p[2]],m->p0[t->p[0]]),
					  sub3D(m->p0[t->p[3]],m->p0[t->p[0]])))/6.0;
	t->young=E;
	t->poisson=nu;
	for(i=0;i<4;i++)
		m->m[t->p[i]]+=V/4.0*rho;
}
void model_newFromMeshFile(Model *m, char *path, double E, double nu, double rho, double thickness)
/*
 Creates a new model from a surface mesh. The mesh is intruded along
 its normal to produce a tetrahedral surface of fixed thickness.
 Furthermore, memory is allocated for all the variables needed for
 the model, and the resting length of radial fibres is configured.
 */
{
	FILE		*f;
	int			*neighb,nverts,ntris,i;
	double3D	*verts,*n,n1;
	int3D		*tris;
	char		str[256];
	float		avrgEdgeLength=0;
	
	// open surface file
	f=fopen(path,"r");
	fscanf(f," %i %i ", &nverts, &ntris);
	verts=(double3D*)calloc(nverts,sizeof(double3D));
	tris=(int3D*)calloc(ntris,sizeof(int3D));
	for(i=0;i<nverts;i++)
	{
		fgets(str,256,f);
		//printf("TEST: %i. %s",i,str);
		sscanf(str," %lf %lf %lf ",&(verts[i].x),&(verts[i].y),&(verts[i].z));
	}
	for(i=0;i<ntris;i++)
	{
		fgets(str,256,f);
		sscanf(str," %i %i %i ",&tris[i].a,&tris[i].b,&tris[i].c);
	}
	fclose(f);
	
	// compute average edge length
	avrgEdgeLength=0;
	for(i=0;i<ntris;i++)
	{
		avrgEdgeLength+=nor3D(sub3D(verts[tris[i].a],verts[tris[i].b]));
		avrgEdgeLength+=nor3D(sub3D(verts[tris[i].b],verts[tris[i].c]));
		avrgEdgeLength+=nor3D(sub3D(verts[tris[i].c],verts[tris[i].a]));
	}
	avrgEdgeLength/=3.0*ntris;	
	
	// tmp for(i=0;i<nverts;i++) verts[i]=sca3D(verts[i],5);

	
	float max[3]={0,0,0},min[3]={0,0,0};
	for(i=0;i<nverts;i++)
	{
		if(verts[i].x>max[0]) max[0]=verts[i].x;
		if(verts[i].x<min[0]) min[0]=verts[i].x;
		if(verts[i].y>max[1]) max[1]=verts[i].y;
		if(verts[i].y<min[1]) min[1]=verts[i].y;
		if(verts[i].z>max[2]) max[2]=verts[i].z;
		if(verts[i].z<min[2]) min[2]=verts[i].z;
	}
	printf("x: %f,%f\ny: %f,%f\nz: %f,%f\n",min[0],max[0],min[1],max[1],min[2],max[2]);
	
	// compute surface normal
	n=(double3D*)calloc(nverts,sizeof(double3D));
	neighb=(int*)calloc(nverts,sizeof(int));
	for(i=0;i<ntris;i++)
	{
		n1=cro3D(sub3D(verts[tris[i].b],verts[tris[i].a]),sub3D(verts[tris[i].c],verts[tris[i].a]));
		n[tris[i].a]=add3D(n[tris[i].a],n1);
		n[tris[i].b]=add3D(n[tris[i].b],n1);
		n[tris[i].c]=add3D(n[tris[i].c],n1);
		neighb[tris[i].a]++;
		neighb[tris[i].b]++;
		neighb[tris[i].c]++;
	}
	for(i=0;i<nverts;i++)
		n[i]=sca3D(n[i],1/(double)neighb[i]);
	free(neighb);
	
	// allocate memory
	m->np=nverts*2;
	m->p0=(double3D*)calloc(m->np,sizeof(double3D));
	m->p=(double3D*)calloc(m->np,sizeof(double3D));
	m->v=(double3D*)calloc(m->np,sizeof(double3D));
	m->m=(double*)calloc(m->np,sizeof(double));
	m->f0=(double3D*)calloc(m->np,sizeof(double3D));
	m->fext=(double3D*)calloc(m->np,sizeof(double3D));
	m->ngb=(Neighbours*)calloc(m->np,sizeof(Neighbours));
	m->fibre_length=(double*)calloc(m->np,sizeof(double));
	m->b=(double3D*)calloc(m->np,sizeof(double3D));
	m->Ap=(double3D*)calloc(m->np,sizeof(double3D));
	m->R=(double3D*)calloc(m->np,sizeof(double3D));
	m->P=(double3D*)calloc(m->np,sizeof(double3D));
	m->nhash=9973; // a prime number close to 10k
	m->hash=(HashCell*)calloc(m->nhash,sizeof(HashCell));
	m->h=avrgEdgeLength;
	
	// make tetrahedral mesh from triangle mesh by extruding
	for(i=0;i<nverts;i++)
	{
		m->p[i*2+0]=m->p0[i*2+0]=verts[i];
		m->p[i*2+1]=m->p0[i*2+1]=sub3D(verts[i],sca3D(n[i],thickness/nor3D(n[i])));
	}
	free(n);
	m->nt=ntris*3;
	m->t=(Tetra*)calloc(m->nt,sizeof(Tetra));
	for(i=0;i<ntris;i++)
	{
		model_setTetra(m,3*i+0,2*tris[i].b+1,2*tris[i].a+1,2*tris[i].c+0,2*tris[i].c+1,E,nu,rho); // b1, a1, c0, c1
		model_setTetra(m,3*i+1,2*tris[i].c+0,2*tris[i].b+1,2*tris[i].a+0,2*tris[i].a+1,E,nu,rho); // c0, b1, a0, a1
		model_setTetra(m,3*i+2,2*tris[i].a+0,2*tris[i].c+0,2*tris[i].b+0,2*tris[i].b+1,E,nu,rho); // a0, c0, b0, b1
	}
	free(verts);
	free(tris);
	
	// init fibres
	for(i=0;i<m->np;i++)
		m->fibre_length[i]=nor3D(m->p0[i]);
}
void model_stiffness(Model *m, Tetra *t)
/*
 Compute stiffness matrices for each pair of vertices in
 tetrahedron t.
 */
{
	int			i,j;
	double3D		x[3];
	double3D		y[4];
	double		det;
	double		V,E,v;
	double		a,b,c;
	
	for(i=1;i<4;i++)
		x[i]=sub3D(m->p0[t->p[i]],m->p0[t->p[0]]);
	det=detMat(vecs2mat(x[1],x[2],x[3]));
	V=det/6.0   /*TEST*/ *(-1);
	
	y[1]=(double3D){	(x[2].y*x[3].z-x[2].z*x[3].y)/det,
		(x[2].z*x[3].x-x[2].x*x[3].z)/det,
		(x[2].x*x[3].y-x[2].y*x[3].x)/det};
	
	y[2]=(double3D){	(x[1].z*x[3].y-x[1].y*x[3].z)/det,
		(x[1].x*x[3].z-x[1].z*x[3].x)/det,
		(x[1].y*x[3].x-x[1].x*x[3].y)/det};
	
	y[3]=(double3D){	(x[1].y*x[2].z-x[1].z*x[2].y)/det,
		(x[1].z*x[2].x-x[1].x*x[2].z)/det,
		(x[1].x*x[2].y-x[1].y*x[2].x)/det};
	
	y[0]=(double3D){ -y[1].x-y[2].x-y[3].x,
		-y[1].y-y[2].y-y[3].y,
		-y[1].z-y[2].z-y[3].z};
	
	E=t->young;
	v=t->poisson;
	a=V*E*(1-v)/(1+v)/(1-2*v);
	b=V*E*v/(1+v)/(1-2*v);
	c=V*E/(1+v) /* TEST */ /2.0;
	
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		{
			t->K[i][j]=(matrix){	a*y[i].x*y[j].x+c*(y[i].y*y[j].y+y[i].z*y[j].z),
				b*y[i].x*y[j].y+c*(y[i].y*y[j].x),
				b*y[i].x*y[j].z+c*(y[i].z*y[j].x),
				b*y[i].y*y[j].x+c*(y[i].x*y[j].y),
				a*y[i].y*y[j].y+c*(y[i].x*y[j].x+y[i].z*y[j].z),
				b*y[i].y*y[j].z+c*(y[i].z*y[j].y),
				b*y[i].z*y[j].x+c*(y[i].x*y[j].z),
				b*y[i].z*y[j].y+c*(y[i].y*y[j].z),
				a*y[i].z*y[j].z+c*(y[i].y*y[j].y+y[i].x*y[j].x)};
		}
}
void model_addEdge(Model *m, int i, int j, Tetra *t)
/*
 This function is used to configure the topology data of
 the model, and is called by model_topology. It find vertex
 t.p[j] among the neighbours of t.p[i].
 */
{
	int	l;
	int	pi,pj;
	int	ij;
	
	pi=t->p[i];
	pj=t->p[j];
	
	for(l=0;l<m->ngb[pi].n;l++)
		if(m->ngb[pi].p[l]==pj)
			break;
	if(l<m->ngb[pi].n)
		ij=l;
	else
		if(m->ngb[pi].n<MAX_NEIG)
		{
			m->ngb[pi].p[m->ngb[pi].n]=pj;
			ij=m->ngb[pi].n;
			m->ngb[pi].n++;
		}
		else
			printf("ERROR: Neighbour overflow for vertex %i\n",i);
	
	t->ngb[i][j]=ij;
}
void model_topology(Model *m)
/*
 Configures the topology data of the model by constructing
 neighbour lists.
 */
{
	int		i,j,l;
	int		min,max;
	
	for(l=0;l<m->nt;l++)
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				model_addEdge(m,i,j,&(m->t[l]));
	min=max=m->ngb[0].n;
	for(i=1;i<m->np;i++)
	{
		min=(m->ngb[i].n<min)?m->ngb[i].n:min;
		max=(m->ngb[i].n>max)?m->ngb[i].n:max;
	}
	printf(" min #ngb: %i\n",min);
	printf(" max #ngb: %i\n",max);
}
#pragma mark -
#pragma mark [ model: simulate ]
void model_rotation(Model *m)
/*
 Computes the rotation from the resting orientation
 of the tetrahedron to its actual orientation. This
 is used for rotation warping in the computation of
 strain (allowing us to use a linear strain tensor
 to handle large deformations)
 */
{
	int		i;
	Tetra	*t;
	matrix	p,x;
	
	for(i=0;i<m->nt;i++)
	{
		t=&(m->t[i]);
		p=vecs2mat( sub3D(m->p[t->p[1]],m->p[t->p[0]]),
				   sub3D(m->p[t->p[2]],m->p[t->p[0]]),
				   sub3D(m->p[t->p[3]],m->p[t->p[0]]));
		x=vecs2mat( sub3D(m->p0[t->p[1]],m->p0[t->p[0]]),
				   sub3D(m->p0[t->p[2]],m->p0[t->p[0]]),
				   sub3D(m->p0[t->p[3]],m->p0[t->p[0]]));
		t->R=ortMat(pinvxMat(p,x));
	}
}
void model_assemble(Model *m)
/*
 Assembles the individual stiffness matrices into a
 single one (sparsely coded)
 */
{
	int			i,j,k,l;
	int			ij,ji;
	double3D		f0;
	Tetra		*t;
	matrix		K;
	
	// f0 = R*K*x
	for(i=0;i<m->np;i++)
		m->f0[i]=(double3D){0,0,0};
	for(i=0;i<m->nt;i++)
	{
		t=&(m->t[i]);
		for(j=0;j<4;j++)
		{
			f0=(double3D){0,0,0};
			for(k=0;k<4;k++)
				f0=add3D(f0,mulMat(t->K[j][k],m->p0[t->p[k]]));
			m->f0[t->p[j]]=add3D(m->f0[t->p[j]],mulMat(t->R,f0));
		}
	}
	
	// K' = R*K*R'
	for(i=0;i<m->np;i++)
		for(j=0;j<m->ngb[i].n;j++)
			m->ngb[i].K[j]=(matrix){0,0,0,0,0,0,0,0,0};
	for(l=0;l<m->nt;l++)
	{
		t=&(m->t[l]);
		for(i=0;i<4;i++)
			for(j=i;j<4;j++)
			{
				K=rkrtMat(t->R,t->K[i][j]);
				ij=t->ngb[i][j];	// j in i's neighbours
				ji=t->ngb[j][i];	// i in j's neighbours
				m->ngb[t->p[i]].K[ij]=addMat(m->ngb[t->p[i]].K[ij],K);		
				if(i!=j)
					m->ngb[t->p[j]].K[ji]=addMat(m->ngb[t->p[j]].K[ji],trnMat(K));
			}
	}
}	
void model_externalForces(Model *m, double plasticityf, double Efib)
/*
 Compute external forces. Here, the action of radial fibres
 */
{
	int			i;
	
	// fibre elasticity
	for(i=1;i<m->np;i+=2)
		m->fext[i]=sca3D(m->p[i],Efib*(m->fibre_length[i]/nor3D(m->p[i])-1));
	
	// fibre plasticity
	for(i=1;i<m->np;i+=2)
		m->fibre_length[i]+=plasticityf*(nor3D(m->p[i])-m->fibre_length[i]);
}
void  model_configureConjugateGradient(Model *m)
/*
 Writes the equation obtained from the implicit Euler scheme
 used to solve the dynamic mechanical equations in the form
 Ax=b, to be solved using the conjugate gradient method
 */
{
	// configure:
	// A =(M-dt^2*K') = (M-dt^2*R*K*R')
	// b = M*v' + dt*(K'p'-f0+fext) = M*v' + dt*(R*K*R'*p'-R*K*v+fext)
	
	int			i,j;
	double3D		b;
	Neighbours	*ngb;
	
	for(i=0;i<m->np;i++)
	{
		ngb=&(m->ngb[i]);
		
		// A = M - dt^2*K'
		for(j=0;j<ngb->n;j++)
		{
			ngb->A[j]=scaMat(ngb->K[j],-m->dt*m->dt);
			if(ngb->p[j]==i)
				ngb->A[j]=addMat(ngb->A[j],(matrix){m->m[i],0,0, 0,m->m[i],0, 0,0,m->m[i]});
		}
		
		// b=M*v +dt*(K'*p-f0+fext)
		b=(double3D){0,0,0};
		for(j=0;j<ngb->n;j++)
			b=add3D(b,mulMat(ngb->K[j],m->p[ngb->p[j]]));
		b=sub3D(b,m->f0[i]);
		b=add3D(b,m->fext[i]);
		b=sca3D(b,m->dt);
		m->b[i]=add3D(sca3D(m->v[i],m->m[i]),b);
	}
}
void model_conjugateGradient(Model *m, int maxiter)
/*
 Conjugate gradient method
 */
{
	int		i,j,k;
	double	den;
	double	rold,rnew;
	double3D	sum;
	double	alpha,beta;
	
	
	// 1. initialise p=r=b-A*x
	for(i=0;i<m->np;i++)
	{
		sum=(double3D){0,0,0};
		for(j=0;j<m->ngb[i].n;j++)
			sum=add3D(sum,mulMat(m->ngb[i].A[j],m->v[m->ngb[i].p[j]]));
		m->P[i]=m->R[i]=sub3D(m->b[i],sum);
	}
	
	// 2. compute the residue rold
	rold=0;
	for(i=0;i<m->np;i++)
		rold+=dot3D(m->R[i],m->R[i]);
	
	for(k=0;k<maxiter;k++)
	{
		// 3. compute A*p
		for(i=0;i<m->np;i++)
		{
			sum=(double3D){0,0,0};
			for(j=0;j<m->ngb[i].n;j++)
				sum=add3D(sum,mulMat(m->ngb[i].A[j],m->P[m->ngb[i].p[j]]));
			m->Ap[i]=sum;
		}
		
		// 4. compute alpha=rold/(p'*A*p)
		den=0;
		for(i=0;i<m->np;i++)
			den+=dot3D(m->P[i],m->Ap[i]);
		if(den<1E-10)
			den=1E-10;
		alpha=rold/den;
		
		// 5. update x and r
		for(i=0;i<m->np;i++)
		{
			m->v[i]=add3D(m->v[i],sca3D(m->P[i],alpha));
			m->R[i]=sub3D(m->R[i],sca3D(m->Ap[i],alpha));
		}
		
		// 6. compute the new residue rnew, finish if error is small enough
		rnew=0;
		for(i=0;i<m->np;i++)
			rnew+=dot3D(m->R[i],m->R[i]);
		if(k>1 && rnew<0.001)
			break;
		
		// 7. compute beta=rnew/rold
		beta=rnew/rold;
		rold=rnew;
		
		// 8. update p
		for(i=0;i<m->np;i++)
			m->P[i]=add3D(m->R[i],sca3D(m->P[i],beta));
	}
	printf("%i %lf ",k,rnew);
}
void model_updatePosition(Model *m)
/*
 Use the velocity computed with the conjugate gradient method
 to update the position of the model vertices
 */
{
	int	i;
	for(i=0;i<m->np;i++)
		m->p[i]=add3D(m->p[i],sca3D(m->v[i],m->dt));
}
void model_addToHash(Model *m, unsigned int hash, int vertexIndex, int iter)
/*
 Hash table 'add' function for collision detection
 */
{
	HashCell	*h=&(m->hash[hash]);
	
	if(h->timeStamp)
	{
		// find a free cell, or add a new one
		while(h->timeStamp==iter)
		{
			if(h->next==NULL)
				h->next=(long*)calloc(1,sizeof(HashCell));
			h=(HashCell*)h->next;
		}
	}
	h->timeStamp=iter;
	h->i=vertexIndex;
}
int model_vertexInTetra(Model *m, int vertexIndex, int tetraIndex, double3D *penetration)
{
	double3D	coords;
	double3D	a,b,c,x;
	
	a=sub3D(m->p[m->t[tetraIndex].p[1]],m->p[m->t[tetraIndex].p[0]]);
	b=sub3D(m->p[m->t[tetraIndex].p[2]],m->p[m->t[tetraIndex].p[0]]);
	c=sub3D(m->p[m->t[tetraIndex].p[3]],m->p[m->t[tetraIndex].p[0]]);
	x=sub3D(m->p[vertexIndex],m->p[m->t[tetraIndex].p[0]]);
	
	coords=invmxp(trnMat(vecs2mat(a, b, c)),x);
	if(coords.x>=0 && coords.y>=0 && coords.z>=0 &&
	   coords.x+coords.y+coords.z<=1 )
		return 1;
	return 0;
}
int model_collision(Model *m,int iter)
/*
 Collision detection algorithm
 */
{
	int				i,j,k,x,y,z;
	unsigned int	hash;
	double3D		min,max,penetration;
	Tetra			*t;
	HashCell		*h;
	int				isInTetrahedron;
	
	// 1st pass: assign vertices to hash
	for(i=0;i<m->np/2;i++)		// np/2 because only the vertices of the external surface are used
	{
		hash=(((unsigned int)floor(m->p[2*i].x/m->h)*92837111)^((unsigned int)floor(m->p[2*i].y/m->h)*689287499)^((unsigned int)floor(m->p[2*i].z/m->h)*283923481))%m->nhash;
		model_addToHash(m,hash,2*i,iter);
	}
	
	// 2nd pass: detect collision with tetrahedra bounding box
	for(i=0;i<m->nt;i++)
	{
		t=&(m->t[i]);
		
		// bounding box for t
		min=max=m->p[t->p[0]];
		for(j=0;j<4;j++)
		{
			if(m->p[t->p[j]].x<min.x) min.x=m->p[t->p[j]].x;
			if(m->p[t->p[j]].y<min.y) min.y=m->p[t->p[j]].y;
			if(m->p[t->p[j]].z<min.z) min.z=m->p[t->p[j]].z;
			
			if(m->p[t->p[j]].x>max.x) max.x=m->p[t->p[j]].x;
			if(m->p[t->p[j]].y>max.y) max.y=m->p[t->p[j]].y;
			if(m->p[t->p[j]].z>max.z) max.z=m->p[t->p[j]].z;
		}
		
		// scan bounding box
		for(x=floor(min.x/m->h);x<=floor(max.x/m->h);x++)
			for(y=floor(min.y/m->h);y<=floor(max.y/m->h);y++)
				for(z=floor(min.z/m->h);z<=floor(max.z/m->h);z++)
				{
					hash=((unsigned int)(x*92837111)^(unsigned int)(y*689287499)^(unsigned int)(z*283923481))%m->nhash;
					
					// check for collision
					h=&(m->hash[hash]);
					while(h->timeStamp)
					{
						if(h->timeStamp==iter)
						{
							// check if the vertex belongs to tetrahedron
							isInTetrahedron=0;
							for(k=0;k<4;k++)
							{
								if(m->t[i].p[k]==h->i)
								{
									isInTetrahedron=1;
									break;
								}
							}
							if(isInTetrahedron==0)
							{
								// 3rd pass: check if vertex is inside tetrahedron
								if(model_vertexInTetra(m,h->i,i,&penetration))
									return 1;
							}
						}
						if(h->next)
							h=(HashCell*)h->next;
						else
							break;
					}
				}
	}
	
	return 0;
}
#pragma mark -
#pragma mark [ util ]
void model_save(Model *m, char *path, int surf)
{
	// save model cortical layer
	
	FILE	*f=fopen(path,"w");
	int		i;
	
	if(surf==1)
	{
		printf("Saving external surface\n");
		fprintf(f,"%i %i\n",m->np/2,m->nt/3);
		for(i=0;i<m->np/2;i++) fprintf(f,"%g %g %g\n",m->p[2*i].x,m->p[2*i].y,m->p[2*i].z);
		for(i=0;i<m->nt/3;i++) fprintf(f,"%i %i %i\n",m->t[3*i].p[1]/2,m->t[3*i].p[0]/2,m->t[3*i].p[2]/2);
	}
	else
		if(surf==2)
		{
			printf("Saving internal surface\n");
			fprintf(f,"%i %i\n",m->np/2,m->nt/3);
			for(i=0;i<m->np/2;i++) fprintf(f,"%g %g %g\n",m->p[2*i+1].x,m->p[2*i+1].y,m->p[2*i+1].z);
			for(i=0;i<m->nt/3;i++) fprintf(f,"%i %i %i\n",m->t[3*i].p[1]/2,m->t[3*i].p[0]/2,m->t[3*i].p[2]/2);
		}
		else
			if(surf==3)
			{
				printf("Saving external and internal surfaces\n");
				fprintf(f,"%i %i\n",m->np,2*m->nt/3);
				for(i=0;i<m->np/2;i++) fprintf(f,"%g %g %g\n",m->p[2*i].x,m->p[2*i].y,m->p[2*i].z);
				for(i=0;i<m->np/2;i++) fprintf(f,"%g %g %g\n",m->p[2*i+1].x,m->p[2*i+1].y,m->p[2*i+1].z);
				for(i=0;i<m->nt/3;i++) fprintf(f,"%i %i %i\n",m->t[3*i].p[1]/2,m->t[3*i].p[0]/2,m->t[3*i].p[2]/2);
				for(i=0;i<m->nt/3;i++) fprintf(f,"%i %i %i\n",m->t[3*i].p[1]/2+m->np/2,m->t[3*i].p[0]/2+m->np/2,m->t[3*i].p[2]/2+m->np/2);
			}
	fclose(f);
}	
void model_free(Model *m)
/*
 Free memory allocated for the model
 */
{
	free(m->p0);
	free(m->p);
	free(m->v);
	free(m->m);
	free(m->f0);
	free(m->fext);
	free(m->ngb);
	free(m->fibre_length);
	free(m->b);
	free(m->Ap);
	free(m->R);
	free(m->P);
	free(m->t);
	
}
double3D	*p;
int3D	*t;
int	np; //number of vertices
int	nt; //number of triangles
double *data;

int curvature(Model *m, double	*C)//chenlu
{
    double3D		*tmp,*tmp1;
    int			*n;
    double3D		nn;
    double		absmax;
    int			i;
	
	np=m->np/2;
	nt=m->nt/3;
	p = (double3D*)calloc(np,sizeof(double3D));
	if(p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
	for(i=0;i<np;i++) {p[i].x=m->p[2*i].x;p[i].y=m->p[2*i].y;p[i].z=m->p[2*i].z;}
	
	t = (int3D*)calloc(nt,sizeof(int3D));
	if(t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
	for(i=0;i<nt;i++) {t[i].a=m->t[3*i].p[1]/2;t[i].b=m->t[3*i].p[0]/2;t[i].c=m->t[3*i].p[2]/2;}
	
	
    tmp=(double3D*)calloc(np,sizeof(double3D));
    n=(int*)calloc(np,sizeof(int));
    // compute smoothing direction as the vector to the average of neighbour vertices
    for(i=0;i<nt;i++)
    {
    	tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
    	tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
    	tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
    	n[t[i].a]+=2;
    	n[t[i].b]+=2;
    	n[t[i].c]+=2;
    }
    for(i=0;i<np;i++)
    	tmp[i]=sub3D(sca3D(tmp[i],1/(double)n[i]),p[i]);
	
    tmp1=(double3D*)calloc(np,sizeof(double3D));
    // compute normal direction as the average of neighbour triangle normals
    for(i=0;i<nt;i++)
    {
    	nn=cro3D(sub3D(p[t[i].b],p[t[i].a]),sub3D(p[t[i].c],p[t[i].a]));
    	nn=sca3D(nn,1/nor3D(nn));
    	tmp1[t[i].a]=add3D(tmp1[t[i].a],nn);
    	tmp1[t[i].b]=add3D(tmp1[t[i].b],nn);
    	tmp1[t[i].c]=add3D(tmp1[t[i].c],nn);
    }
    for(i=0;i<np;i++)
    	tmp1[i]=sca3D(tmp1[i],1/(double)n[i]);
    free(n);
    
    for(i=0;i<np;i++)
		C[i]=-dot3D(tmp1[i],tmp[i]);
	free(tmp);
	free(tmp1);
    absmax=-1;
    for(i=0;i<np;i++)
        absmax=(fabs(C[i])>absmax)?fabs(C[i]):absmax;
    for(i=0;i<np;i++)
        C[i]/=absmax;
    
    return 0;
}
void foldLength(void)//chenlu
{
	int		i,j;
	double	length=0,a,x;
	double3D	p0[3];
	
	for(i=0;i<nt;i++)
	{
		j=0;
		if(data[t[i].a]*data[t[i].b]<0)
		{
			a=fabs(data[t[i].a]);
			x=a/(a+fabs(data[t[i].b]));
			p0[j++]=add3D(sca3D(p[t[i].a],1-x),sca3D(p[t[i].b],x));
		}
		if(data[t[i].b]*data[t[i].c]<0)
		{
			a=fabs(data[t[i].b]);
			x=a/(a+fabs(data[t[i].c]));
			p0[j++]=add3D(sca3D(p[t[i].b],1-x),sca3D(p[t[i].c],x));
		}
		if(data[t[i].c]*data[t[i].a]<0)
		{
			a=fabs(data[t[i].c]);
			x=a/(a+fabs(data[t[i].a]));
			p0[j++]=add3D(sca3D(p[t[i].c],1-x),sca3D(p[t[i].a],x));
		}
		if(j==2)
			length+=nor3D(sub3D(p0[0],p0[1]));
	}
	printf("%g\n",length/2.0);
}

#pragma mark -
int main(int argc, char *argv[])
/*
 Make a model from a mesh, initialise it,
 simulate it, save the result, free it
 
 arg[1] is the mesh that gives the shape of the model surface
 arg[2] is the files where the deformed mesh will be saved
 arg[3]	E: Young modulus
 arg[4]	nu: Poisson ratio
 arg[5]	rho: density of the material (mass/volume)
 arg[6]	thickness of the cortical layer
 arg[7]	number of iterations
 arg[8] surfaces to save: 1=external surface (default), 2=internal surface, 3=both (default)
 
 */
{
	Model	m;
	int		i,j;
	char	*input,*output;
	int		niter=100;
	float	E=400000;
	float	nu=0.33;
	float	rho=1000;
	float	thickness=0.4;
	int		surf=1;
	int		step=0;
	float	tgrowth=0.1;
	float	plasticityf=0.00001;
	float	Efib=1400;
	float	plasticity=0.006125;
	char	str[1024];
	double	r;
	double	initialVolume,actualVolume,targetVolume;
	
	
	printf("#arguments=%i\n",argc);
	
	// read arguments
	i=1;
	while(i<argc)
	{
		if(strcmp(argv[i],"-i")==0)
			input=argv[++i];
		else
		if(strcmp(argv[i],"-o")==0)
			output=argv[++i];
		else
		if(strcmp(argv[i],"-E")==0)
			E=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-nu")==0)
			nu=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-rho")==0)
			rho=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-thickness")==0)
			thickness=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-plasticity")==0)
			plasticity=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-plasticityf")==0)
			plasticityf=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-Efib")==0)
			Efib=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-tgrowth")==0)
			tgrowth=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-niter")==0)
			niter=atof(argv[++i]);
		else
		if(strcmp(argv[i],"-surf")==0)
			surf=atoi(argv[++i]);
		else
		if(strcmp(argv[i],"-step")==0)
			step=atoi(argv[++i]);
		i++;
	}
	
	// load mesh
	printf("load mesh\n");
	model_newFromMeshFile(&m,input,E,nu,rho,thickness);
	//model_newFromMeshFile(&m,argv[1],400000,0.33,1000,0.1);
	//model_newFromMeshFile(&m,argv[1],4000000,0.33,1000,0.05);
	
	
	// init element stiffness matrices
	printf("init stiffness matrices\n");
	for(i=0;i<m.nt;i++)							// K
		model_stiffness(&m, &(m.t[i]));
	
	// init mesh topology
	printf("init model topology\n");
	model_topology(&m);
	
	
	// growth
	initialVolume=0;
	for(j=0;j<m.nt;j++)
		initialVolume+=detMat(vecs2mat(sub3D(m.p0[m.t[j].p[1]],m.p0[m.t[j].p[0]]),
									  sub3D(m.p0[m.t[j].p[2]],m.p0[m.t[j].p[0]]),
									  sub3D(m.p0[m.t[j].p[3]],m.p0[m.t[j].p[0]])))/6.0;
	printf("Initial volume: %lf\n", initialVolume);
	targetVolume=4*initialVolume;
	printf("Target volume: %lf\n", targetVolume);
	
	
	// simulation loop
	printf("start simulation loop\n");
	m.dt=0.1;
	for(i=0;i<niter;i++)
	{
		printf("%i ",i);
		
		// cortical layer growth
		actualVolume=0;
		for(j=0;j<m.nt;j++)
		{
			actualVolume+=detMat(vecs2mat(sub3D(m.p0[m.t[j].p[1]],m.p0[m.t[j].p[0]]),
										 sub3D(m.p0[m.t[j].p[2]],m.p0[m.t[j].p[0]]),
										 sub3D(m.p0[m.t[j].p[3]],m.p0[m.t[j].p[0]])))/6.0;
		}
		
		
		r=pow(tgrowth*(actualVolume/initialVolume-0.99)*pow(1-actualVolume/targetVolume,2)+1,1/3.0);
		printf("%lf %lf ",actualVolume,r);
		for(j=0;j<m.np;j++)
			m.p0[j]=sca3D(m.p0[j],r);
		
		model_rotation(&m);						// compute R
		model_assemble(&m);						// compute f0=R*K, K'=R*K*R'
		model_externalForces(&m,plasticityf,Efib);			// compute fext
		model_configureConjugateGradient(&m);	// compute A=M-dt^2*K', b=M*v-dt*(K'*p-f0+fext)
		model_conjugateGradient(&m,10);			// solve A*v=b for the vertex velocities v
		model_updatePosition(&m);				// update vertex position based on velocities
		
		// cortical layer plasticity
		for(j=0;j<m.np;j++)
			m.p0[j]=add3D(m.p0[j],sca3D(sub3D(m.p[j],m.p0[j]),plasticity));
		
		// TEST: make rest configuration equals to the actual
		//for(j=0;j<m.np;j++) m.p0[j]=m.p[j];
		
		// TEST: recompute stiffness matrices
		//for(j=0;j<m.nt;j++) model_stiffness(&m, &(m.t[j]));
		
		if(step>0 && i%step==0)
		{
			sprintf(str,"%s.%i.txt",output,i);
			model_save(&m,str,surf);
		}
		
		data=NULL;//chenlu
		if(data==NULL)
			data=(double*)calloc(m.np/2,sizeof(double));
		curvature(&m,data);
		foldLength();
		if(data)//chenlu
			free(data);	

		if(model_collision(&m,i+1))
		{
			printf("Collision detected.\n");
			sprintf(str,"%s.%i.txt",output,i);
			model_save(&m,str,surf);	
			return 0;
		}		
	}

	actualVolume=0;
	for(j=0;j<m.nt;j++)
		actualVolume+=detMat(vecs2mat(sub3D(m.p0[m.t[j].p[1]],m.p0[m.t[j].p[0]]),
									  sub3D(m.p0[m.t[j].p[2]],m.p0[m.t[j].p[0]]),
									  sub3D(m.p0[m.t[j].p[3]],m.p0[m.t[j].p[0]])))/6.0;
	printf("Final volume: %lf\n",actualVolume);
	
	sprintf(str,"%s.%i.txt",output,i);
	model_save(&m,str,surf);
	
	model_free(&m);
	
	
	return 0;
}

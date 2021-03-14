#ifndef AMREX_EB_TRIGEOMOPS_K_H_
#define AMREX_EB_TRIGEOMOPS_K_H_
#include <AMReX_Config.H>
#include <AMReX.H>

namespace amrex
{
    namespace tri_geom_ops
    {
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real Distance2(Real P1[3],Real P2[3])
        {
            return( (P1[0]-P2[0])*(P1[0]-P2[0]) +
                    (P1[1]-P2[1])*(P1[1]-P2[1]) +
                    (P1[2]-P2[2])*(P1[2]-P2[2]) );
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real DotProd(Real v1[3],Real v2[3])
        {
            return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real side_op(Real L1[6],Real L2[6])
        {
            return(    L1[0]*L2[4]
                    +  L1[1]*L2[5]
                    +  L1[2]*L2[3]
                    +  L1[3]*L2[2]
                    +  L1[4]*L2[0]
                    +  L1[5]*L2[1] );
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void getvec(Real P1[3],Real P2[3],Real v[3])
        {
            v[0]=P2[0]-P1[0];
            v[1]=P2[1]-P1[1];
            v[2]=P2[2]-P1[2];
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void getunitvec(Real v[3],Real vu[3])
        {
            Real vmag;
            vmag=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            vu[0]=v[0]/vmag;
            vu[1]=v[1]/vmag;
            vu[2]=v[2]/vmag;
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CrossProd(Real v1[3],Real v2[3],Real v[3])
        {
            v[0]=v1[1]*v2[2]-v1[2]*v2[1];
            v[1]=v1[2]*v2[0]-v1[0]*v2[2];
            v[2]=v1[0]*v2[1]-v1[1]*v2[0];
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_plucker_coords(Real v1[3],Real v2[3],Real L[6])
        {
            L[0] = v1[0]*v2[1] - v1[1]*v2[0];
            L[1] = v1[0]*v2[2] - v1[2]*v2[0];
            L[2] = v1[0]       - v2[0];
            L[3] = v1[1]*v2[2] - v1[2]*v2[1];
            L[4] = v1[2]       - v2[2];
            L[5] = v2[1]       - v1[1];
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void side_op3(Real v1[3],Real v2[3],
                Real t1[3],Real t2[3],Real t3[3],
                Real &S1, Real &S2, Real &S3)
        {

            Real L[6],e1[6],e2[6],e3[6];

            get_plucker_coords(v1,v2,L);
            get_plucker_coords(t1,t2,e1);
            get_plucker_coords(t2,t3,e2);
            get_plucker_coords(t3,t1,e3);

            S1=side_op(L,e1);
            S2=side_op(L,e2);
            S3=side_op(L,e3);
        }
        //================================================================================
        //get normal of triangle pointing at a test-point
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void tri_n(Real P1[3],Real P2[3],Real P3[3],
                Real testp[3],Real n[3])
        {
            Real v1[3],v2[3],magn;
            Real centr[3],c_tp_vec[3];

            getvec(P1,P2,v1);
            getvec(P1,P3,v2);
            CrossProd(v1,v2,n);


            centr[0]=0.333333*(P1[0]+P2[0]+P3[0]);
            centr[1]=0.333333*(P1[1]+P2[1]+P3[1]);
            centr[2]=0.333333*(P1[2]+P2[2]+P3[2]);

            getvec(centr,testp,c_tp_vec);
            magn=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

            if(DotProd(c_tp_vec,n) < 0.0)
            {
                magn=-magn;
            }

            n[0]=n[0]/magn;
            n[1]=n[1]/magn;
            n[2]=n[2]/magn;
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real triangle_area(Real P1[3],Real P2[3],Real P3[3])
        {
            Real v1[3],v2[3],area[3];

            getvec(P1,P2,v1);
            getvec(P1,P3,v2);
            CrossProd(v1,v2,area);
            return(0.5 * std::sqrt(DotProd(area,area)) );
        }
        //================================================================================
        //this is only useful when v1-v2 segment intersects the triangle
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool find_intersection_point(Real v1[3],Real v2[3],
                Real t1[3], Real t2[3], Real t3[3],Real ip[3],int bisect_iters=20,Real tol=1e-6)
        {
            Real plane_eq_mid,plane_eq1,plane_eq2;

            Real ab[3],ac[3],n[3],magn;
            Real midp[3],p1[3],p2[3];

            getvec(t1,t2,ab);
            getvec(t1,t3,ac);

            CrossProd(ab,ac,n);
            magn=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

            n[0]=n[0]/magn;
            n[1]=n[1]/magn;
            n[2]=n[2]/magn;

            p1[0]=v1[0];
            p1[1]=v1[1];
            p1[2]=v1[2];

            p2[0]=v2[0];
            p2[1]=v2[1];
            p2[2]=v2[2];


            bool all_ok=true;

            for(int i=0;i<bisect_iters;i++)
            {
                midp[0]=0.5*(p1[0]+p2[0]);
                midp[1]=0.5*(p1[1]+p2[1]);
                midp[2]=0.5*(p1[2]+p2[2]);

                plane_eq_mid= (midp[0]-t1[0])*n[0] + (midp[1]-t1[1])*n[1] + (midp[2]-t1[2])*n[2];
                plane_eq1   = (p1[0]  -t1[0])*n[0] + (p1[1]  -t1[1])*n[1] + (p1[2]  -t1[2])*n[2];
                plane_eq2   = (p2[0]  -t1[0])*n[0] + (p2[1]  -t1[1])*n[1] + (p2[2]  -t1[2])*n[2];

                //Print()<<"midp:"<<midp[0]<<"\t"<<midp[1]<<"\t"<<midp[2]<<"\t"<<plane_eq_mid<<"\n";

                if(Math::abs(plane_eq_mid) < tol)
                {
                    break;
                }

                if(plane_eq_mid*plane_eq1 < 0.0)
                {
                    p2[0]=midp[0];
                    p2[1]=midp[1];
                    p2[2]=midp[2];
                }
                else if(plane_eq_mid*plane_eq2 < 0.0)
                {
                    p1[0]=midp[0];
                    p1[1]=midp[1];
                    p1[2]=midp[2];
                }
                else //plane_eq_mid is 0
                    //or error: p1,midp and p2 are on the same side
                    //which is not what this function is meant for
                {
                    if(plane_eq_mid*plane_eq1 > 0.0 && plane_eq_mid*plane_eq2 > 0.0)
                    {
                        all_ok=false;
                    }
                    break;
                }
            }

            ip[0]=midp[0];
            ip[1]=midp[1];
            ip[2]=midp[2];

            return(all_ok);
        }
        //================================================================================
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int lineseg_tri_intersect(Real v1[3],Real v2[3],
                Real t1[3],Real t2[3],Real t3[3])
        {
            //see plucker coordinates based method
            //https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/Line_Triangle.html

            Real S1,S2,S3;
            Real tri_area,area1,area2;
            Real L2[6],L3[6],L4[6],ls_s1,ls_s2;

            side_op3(v1,v2,t1,t2,t3,S1,S2,S3);

            //we are assuming there are no intersections initially
            int no_intersections=1;

            Real eps = std::numeric_limits<Real>::epsilon();

            //coplanar (S1,S2,S3 = 0)
            if(Math::abs(S1) < eps && Math::abs(S2) < eps && Math::abs(S3) < eps)
            {
                //Print()<<"line segment and triangle are in the same plane\t"<<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
                tri_area=triangle_area(t1,t2,t3);

                /*if(tri_area == 0)
                  {
                  amrex::Abort("problem with triangle\n");
                  }*/
                area1=(triangle_area(t1,t2,v1)+triangle_area(t2,t3,v1)+triangle_area(t3,t1,v1));
                area2=(triangle_area(t1,t2,v2)+triangle_area(t2,t3,v2)+triangle_area(t3,t1,v2));

                if( Math::abs(area1-tri_area)>eps || Math::abs(area2-tri_area)>eps )
                {
                    no_intersections = 0;
                }
            }
            //proper and edge intersection
            else if( (S1 < 0.0 && S2 < 0.0 && S3 < 0.0) ||
                    (S1 > 0.0 && S2 > 0.0 && S3 > 0.0) ||
                    (Math::abs(S1) < eps && S2*S3 > 0.0) ||     //S1=0
                    (Math::abs(S2) < eps && S3*S1 > 0.0) ||     //S2=0
                    (Math::abs(S3) < eps && S1*S2 > 0.0) )      //S3=0
            {

                get_plucker_coords(v1,t1,L2);
                get_plucker_coords(t1,v2,L3);
                get_plucker_coords(t2,t3,L4);

                /*if(Math::abs(S1*S2*S3) < eps)
                  {
                  Print()<<"edge intersection S1,S2,S3:"
                  <<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
                  }*/

                ls_s1 = side_op(L4,L3);
                ls_s2 = side_op(L4,L2);

                if(ls_s1*ls_s2 > 0.0)
                {
                    no_intersections = 0;
                }
            }
            //vertex intersection
            else if( (Math::abs(S1) < eps && Math::abs(S2) < eps) ||  //S1,S2=0
                    (Math::abs(S2) < eps && Math::abs(S3) < eps) )   //S2,S3=0
            {

                //Print()<<"vertex intersection type 1\n";
                //don't chose vertex 2 or 3
                get_plucker_coords(v2,t1,L2);
                get_plucker_coords(t1,v1,L3);
                get_plucker_coords(t2,t3,L4);

                ls_s1=side_op(L4,L3);
                ls_s2=side_op(L4,L2);

                if(ls_s1*ls_s2 > 0)
                {
                    no_intersections = 0;
                }
            }
            else if(Math::abs(S3) < eps && Math::abs(S1) < eps) //S3,S1=0
            {

                //Print()<<"vertex intersection type 2\n";
                //don't chose vertex 1
                get_plucker_coords(v2,t2,L2);
                get_plucker_coords(t2,v1,L3);
                get_plucker_coords(t3,t1,L4);

                ls_s1=side_op(L4,L3);
                ls_s2=side_op(L4,L2);

                if(ls_s1*ls_s2 > 0)
                {
                    no_intersections=0;
                }
            }

            return(no_intersections);

        }
        //================================================================================
    }
}
#endif

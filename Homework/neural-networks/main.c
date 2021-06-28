#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

typedef struct {double(*f)(double); double(*f_diff)(double); double(*f_int)(double); gsl_vector* params; } ann;

static int N_MAX;

// Minimization functions
void gradient(double f(gsl_vector* x),gsl_vector* x, gsl_vector* gradient,
              gsl_vector* Dx){
    int n = x -> size;
    gsl_vector_memcpy(Dx,x);
    for (int j = 0; j < n; j++) {
        double xj = gsl_vector_get(x, j);
        gsl_vector_set(Dx, j, xj+sqrt(DBL_EPSILON));
        double J_j = (f(Dx)-f(x))/sqrt(DBL_EPSILON);
        gsl_vector_set(gradient, j, J_j);
        gsl_vector_memcpy(Dx,x);
    }
}
void qnewton( double f(gsl_vector* x),
              gsl_vector* x, double eps){
    int n = x->size;
    double alpha = 1e-4;
    double prikprod;
    double prikprodsy;
    gsl_matrix* B = gsl_matrix_alloc(n,n);
    gsl_vector* grad = gsl_vector_alloc(n);
    gsl_vector* xs = gsl_vector_alloc(n);
    gsl_vector* grad_vec = gsl_vector_alloc(n);
    gsl_vector* s = gsl_vector_alloc(n);
    gsl_vector* grad_xs=gsl_vector_alloc(n);
    gsl_vector* y=gsl_vector_alloc(n);
    gsl_vector* u=gsl_vector_alloc(n);
    N_MAX = 0;
    gsl_matrix_set_identity(B);
    gradient(f,x,grad,grad_vec);
    while(gsl_blas_dnrm2(grad)>eps && N_MAX<10000){
        N_MAX++;
        double lambda = 1;
        gsl_blas_dgemv(CblasNoTrans,-1.,B,grad,0.,s);
        for(int i=0; i<n; i++){
            gsl_vector_set(xs,i, gsl_vector_get(s,i)+gsl_vector_get(x,i));
        }
        gsl_vector_memcpy(xs,x);
        gsl_vector_add(xs,s);
        gsl_blas_ddot(s,grad,&prikprod);
        double fx = f(x);
        while(f(xs)>fx+alpha*prikprod){
            if(lambda<sqrt(DBL_EPSILON)){
                gsl_matrix_set_identity(B);
                break;
            }
            lambda/=2;
            gsl_vector_scale(s,0.5);
            gsl_vector_memcpy(xs,x);
            gsl_vector_add(xs,s);
            gsl_blas_ddot(s,grad,&prikprod);
        }
        gradient(f,xs,grad_xs,grad_vec);
        gsl_vector_memcpy(y,grad_xs);
        gsl_blas_daxpy(-1,grad,y);
        gsl_vector_memcpy(u,s);
        gsl_blas_dgemv(CblasNoTrans,-1.,B,y,1.,u);
        gsl_blas_ddot(s,y,&prikprodsy);
        if(fabs(prikprodsy)>1e-12){
            gsl_blas_dger(1./prikprodsy,u,s,B);
        }
        gsl_vector_memcpy(x,xs);
        gsl_vector_memcpy(grad ,grad_xs);
    }
    gsl_vector_free(grad);
    gsl_matrix_free(B);
    gsl_vector_free(grad_vec);
    gsl_vector_free(xs);
    gsl_vector_free(y);
    gsl_vector_free(grad_xs);
    gsl_vector_free(u);
    gsl_vector_free(s);
}

// Amoeba Functions
void simplex_update(gsl_matrix* simplex, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    int n = centroid -> size;
    *hi = 0;
    *lo = 0;
    double highest = gsl_vector_get(F_val,*hi);
    double lowest = gsl_vector_get(F_val,*lo);
    for(int i = 1; i<n+1;i++){
        double fi = gsl_vector_get(F_val,i);
        if(fi>highest){
            *hi = i;
            highest = fi;
        }
        else if(fi<lowest){
            *lo = i;
            lowest = fi;
        }
    }

    for(int i =0; i<n; i++){
        double sum = 0;
        for(int j =0; j<n+1;j++){
            if(j != *hi){
                sum+=gsl_matrix_get(simplex,i,j)/n;
            }
            gsl_vector_set(centroid,i,sum);
        }
    }
}

void simplex_init(double cost(ann* network, gsl_vector* xs, gsl_vector* ys), gsl_matrix* simplex, ann* init,
               gsl_vector* xs, gsl_vector* ys, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    int n = centroid -> size;
    for(int i = 0; i<n+1; i++){
        gsl_vector_view x = gsl_matrix_column(simplex,i);
        init->params = &x.vector;
        gsl_vector_set(F_val,i,cost(init,xs,ys));
    }
    simplex_update(simplex,F_val,centroid,hi,lo);
}

double size(gsl_matrix* simplex, int lo){
    double s=0.;
    double distance;
    int n = simplex -> size1;
    gsl_vector* dist = gsl_vector_alloc(n);
    gsl_vector_view lowest = gsl_matrix_column(simplex,lo);
    for(int i =0; i<n+1;i++){
        gsl_vector_view xi = gsl_matrix_column(simplex,i);
        gsl_vector_memcpy(dist,&lowest.vector);
        gsl_vector_sub(dist,&xi.vector);
        distance = gsl_blas_dnrm2(dist);
        if(distance>s){
            s=distance;
        }
    }
    gsl_vector_free(dist);
    return s;
}

void reflection(gsl_vector* highest, gsl_vector* centroid, gsl_vector* reflected){
    gsl_vector_memcpy(reflected,highest);
    gsl_vector_scale(reflected,-1.);
    gsl_blas_daxpy(2.,centroid,reflected);
}
void expansion(gsl_vector* highest, gsl_vector* centroid, gsl_vector* expanded){
    gsl_vector_memcpy(expanded,highest);
    gsl_vector_scale(expanded,-2.);
    gsl_blas_daxpy(3.,centroid,expanded);
}
void contraction(gsl_vector* highest, gsl_vector* centroid, gsl_vector* contracted){
    gsl_vector_memcpy(contracted,highest);
    gsl_vector_scale(contracted,0.5);
    gsl_blas_daxpy(0.5,centroid,contracted);
}

void reduction(gsl_matrix* simplex, int lo){
    int n = simplex -> size1;
    for(int i=0; i<n+1;i++){
        if(i!=lo){
            gsl_vector_view lowest = gsl_matrix_column(simplex,lo);
            gsl_vector_view i_vec = gsl_matrix_column(simplex,i);
            gsl_vector_add(&i_vec.vector,&lowest.vector);
            gsl_vector_scale(&i_vec.vector,0.5);
        }
    }
}

void amoeba( double cost(ann* network, gsl_vector* xs, gsl_vector* ys),
             ann* network, gsl_vector* xs, gsl_vector* ys, double eps){
    int n = (network->params)->size;
    gsl_vector* step = gsl_vector_calloc(n);
    gsl_vector_memcpy(step,(network->params));
    for(int i=0; i<n; i++){
        double stepi = gsl_vector_get(step,i);
        gsl_vector_set(step,i,0.7*stepi);
    }
    N_MAX = 0;
    gsl_matrix* simplex = gsl_matrix_alloc(n,n+1);
    gsl_vector* centroid = gsl_vector_alloc(n);

    gsl_vector* init_vec = gsl_vector_alloc(n);
    ann* init = malloc(sizeof(ann));
    init->f = network->f;
    init->params = init_vec;

    ann* ann_p1 = malloc(sizeof(ann));
    gsl_vector* p1 = gsl_vector_alloc(n);
    ann_p1->f = (network->f);
    ann_p1->params = p1;

    ann* ann_p2 = malloc(sizeof(ann));
    gsl_vector* p2 = gsl_vector_alloc(n);
    ann_p2->f = (network->f);
    ann_p2 -> params = p2;

    gsl_vector* F_val = gsl_vector_alloc(n+1);
    int hi, lo;

    for(int i = 0; i<n+1; i++){
        gsl_matrix_set_col(simplex,i,(network->params));
    }
    for(int i = 0; i<n; i++){
        double ii = gsl_matrix_get(simplex,i,i);
        double stepi = gsl_vector_get(step,i);
        gsl_matrix_set(simplex,i,i, ii+stepi);
    }

    simplex_init(cost,simplex,init, xs, ys,F_val,centroid,&hi,&lo);

    while(size(simplex,lo)>eps && N_MAX<1e7){
        N_MAX++;
        simplex_update(simplex,F_val,centroid,&hi,&lo);
        gsl_vector_view highest = gsl_matrix_column(simplex,hi);

    reflection(&highest.vector, centroid, p1);
    double f_re = cost(ann_p1,xs,ys);
        
        if(f_re< gsl_vector_get(F_val,lo)){
            expansion(&highest.vector,centroid,p2);
            double f_ex = cost(ann_p2,xs,ys);
            if(f_ex<f_re){
                gsl_vector_memcpy(&highest.vector,p2);
                gsl_vector_set(F_val,hi,f_ex);
            }
            else{
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }
        }
        else{
            if(f_re<gsl_vector_get(F_val,hi)){
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }

            else {
                contraction(&highest.vector, centroid, p1);
                double f_co = cost(ann_p1,xs,ys);
                if (f_co < gsl_vector_get(F_val, hi)) {
                    gsl_vector_memcpy(&highest.vector, p1);
                    gsl_vector_set(F_val, hi, f_co);
                }
                else {
                    reduction(simplex, lo);
                    simplex_init(cost,simplex,init, xs, ys,F_val,centroid,&hi,&lo);
                }
            }

        }
    }
    gsl_vector_view final_low = gsl_matrix_column(simplex,lo);
    gsl_vector_memcpy((network->params),&final_low.vector);
    gsl_vector_free(centroid);
    gsl_vector_free(F_val);
    gsl_vector_free(p1);
    gsl_vector_free(p2);
    free(ann_p1);
    free(ann_p2);
    free(init);
    gsl_matrix_free(simplex);
}

//Neural Network Functions
ann*   ann_alloc   (int n,double(*f)(double),double(*f_diff)(double),double(* f_int)(double )){
    ann* network = malloc(sizeof(ann));
    gsl_vector* params = gsl_vector_alloc(3*n);
    network->params = params;
    network->f =f;
    network->f_diff = f_diff;
    network -> f_int = f_int;
    return network;
}

void   ann_free    (ann* network){
    gsl_vector_free(network->params);
}

double ann_response(ann* network,double x){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=((network->f)((x-ai)/bi))*wi;
    }
    return Fp;
}

double ann_diff(ann* network, double x){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=((network->f_diff)((x-ai)/bi))*wi/bi;
    }
    return Fp;
}

double ann_int(ann* network, double a, double b){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=(((network->f_int)((b-ai)/bi))*wi*bi)-(((network->f_int)((a-ai)/bi))*wi*bi);
    }
    return Fp;
}

double cost_function(ann* network, gsl_vector* xs, gsl_vector* ys){
    int N = xs->size;
    double cost = 0;
    for(int i=0; i<N; i++){
        double xi = gsl_vector_get(xs,i);
        double yi = gsl_vector_get(ys,i);
        cost += pow((ann_response(network,xi)-yi),2);
    }
    return cost/N;
}

void   ann_train   (ann* network,gsl_vector* xs,gsl_vector* ys){
    amoeba(cost_function,network,xs,ys,1e-7);
}
//Activation functions
double gaussWave(double x){
    return x*exp(-x*x);
}

double gaussWave_diff(double x) {
    return exp(-x * x) - 2 * x * x * exp(-x * x);
}

double gaussWave_int(double x){
    return -0.5*exp(-x*x);
}


// Functions to interpolate
double gauss(double x){
    return exp(-x*x);
}
double gauss_int(double a, double b){
    return 0.5*sqrt(M_PI)*(erf(b)-erf(a));
}
double gauss_diff(double x) {
    return -2 * x * exp(-x * x);
}

int main(){
    printf("-------------- Execise A+B ------------------\n");
    gsl_vector* xs = gsl_vector_alloc(20);
    gsl_vector* ys = gsl_vector_alloc(20);
    printf("#Interpolation data\n");
    for(int i =-10; i<10;i++){
            double xi = (double) i/2;
            double yi = gauss(xi);
            gsl_vector_set(xs,i+10,xi);
            gsl_vector_set(ys,i+10,yi);
        printf("%10g %10g\n", xi, yi);
        }
    ann* interp = ann_alloc(3,gaussWave, gaussWave_diff, gaussWave_int);
    printf("\n\n\n");
    //double a = -5, b= 10;
    int n_neuron_interp = ((interp->params)->size)/3;

    for(int i=0; i<n_neuron_interp;i++){
        //double ai = (b_interp-a_interp)*(i)/(n_neuron_interp-1);
        double bi = 1;
        double wi = 1;

        gsl_vector_set((interp->params),3*i,1);
        gsl_vector_set((interp->params),3*i+1,bi);
        gsl_vector_set((interp->params),3*i+2,wi);
    }

    ann_train(interp,xs,ys);
    printf("# Fitting Data\n");

    for(int i=0; i<400; i++){
        double xi = -2+(double) i/100;
        double interp_i = ann_response(interp,xi);
        double interp_int_i = ann_int(interp,-5,xi);
        double interp_diff_i = ann_diff(interp,xi);
        double exact_diff_i = gauss_diff(xi);
        double exact_int_i = gauss_int(-5,xi);
        double exact_i = gauss(xi);
        printf("%g %g %g %g %g %g %g\n",xi,interp_i,exact_i,interp_diff_i,exact_diff_i,
                interp_int_i,exact_int_i);
    }
    return 0;
}


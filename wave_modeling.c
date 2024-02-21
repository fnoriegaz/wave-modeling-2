#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void spatial_diff(float *field_1, float *field_0, int *model_size, float delta, char axis){
    int jump;
    if(axis=='x'){
        jump = 1;
    }
    else if(axis=='z'){
        jump = model_size[1];
    }
    
    for(int x=1;x<model_size[0]-1;x++){
        for(int z=0;z<model_size[1]-1;z++){
            int idx = x+z*model_size[0];
            field_1[idx] = (field_0[idx+jump]/2.0 - field_0[idx-jump]/2.0)/delta;
        }
    }
}


int main(){

    int model_size[2] = {512,512};
    int source_loc[2] = {512/2,512/2};
    int time_samples = 3000;
    float dt = 1e-3;
    float dx = 12.5;
    float dz = 12.5;
    float freq = 3.0;
    float a = M_PI*M_PI*freq*freq;

    float *wavelet = malloc(sizeof(float)*time_samples);;
    float *c = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *rho = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vx = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vz = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *P_0 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *P_1 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *P_2 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *temp_field = malloc(sizeof(float)*model_size[0]*model_size[1]);


    //creating wavelet
    for(int t=0;t<time_samples;t++){
        float A = (1-2*a*t*t*dt*dt)*exp(-a*t*t*dt*dt);
        wavelet[t] = A;
    }

    //setting rho and c to fixed values
    for(int x=0;x<model_size[0]*model_size[1];x++){
        c[x] = 1540.0;
        rho[x] = 1.0;
    }
    for(int t=0;t<time_samples;t++){
        temp_field = P_2;
        P_2 = P_1;
        P_1 = P_0;
        
        memcpy(temp_field, vx, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(vx, temp_field, model_size, dx, 'x');

        memcpy(temp_field, vz, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(vz, temp_field, model_size, dz, 'z');

        P_0 = temp_field;
        for(int x=1;x<model_size[0]-1;x++){
            for(int z=0;z<model_size[1]-1;z++){
                int idx = x+z*model_size[0];
                P_0[idx] = 2.0 * ( dt*rho[idx]*c[idx]*c[idx]*( vx[idx] + vz[idx] ) + P_2[idx]/2.0 );
            }
        }
        P_0[source_loc[0] + model_size[0]*source_loc[1]] += wavelet[t];

        //TODO: implement vx and vz the same way as P
        memcpy(temp_field, P_0, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(P_0, temp_field, model_size, dx, 'x');
        for(int x=1;x<model_size[0]-1;x++){
            for(int z=0;z<model_size[1]-1;z++){
                int idx = x+z*model_size[0];
                P_0[idx] = 2.0 * ( dt*rho[idx]*c[idx]*c[idx]*( vx[idx] + vz[idx] ) + P_2[idx]/2.0 );
            }
        }

    }

    


    free(c);
    free(rho);
    free(vx);
    free(vz);

    return 0;
}
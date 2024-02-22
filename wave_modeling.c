#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <SDL2/SDL.h>
#include <stdbool.h>

#define SCREEN_WIDTH 768
#define SCREEN_HEIGHT 768

void spatial_diff(float *field_1, float *field_0, int *model_size, float delta, char axis){
    int jump;
    if(axis=='x'){
        jump = 1;
    }
    else if(axis=='z'){
        jump = model_size[1];
    }
    
    for(int x=1;x<model_size[0]-1;x++){
        for(int z=1;z<model_size[1]-1;z++){
            int idx = x+z*model_size[0];
            field_1[idx] = (field_0[idx+jump]/2.0 - field_0[idx-jump]/2.0)/delta;
        }
    }
}

void putpixel(SDL_Surface *surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to set */
    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

    switch(bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16 *)p = pixel;
        break;

    case 3:
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        } else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32 *)p = pixel;
        break;
    }
}

int main(){

    //The window we'll be rendering to
    SDL_Window* window = NULL;
    
    //The surface contained by the window
    SDL_Surface* screenSurface = NULL;

    //Initialize SDL
    if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
    {
        printf( "SDL could not initialize! SDL_Error: %s\n", SDL_GetError() );
        return -1;
    }

    //Create window
    window = SDL_CreateWindow( "SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
    if( window == NULL ){
        printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
        return -1;
    }
    //Get window surface
    screenSurface = SDL_GetWindowSurface( window );
    

    int model_size[2] = {512,512};
    int source_loc[2] = {512/2,512/2};
    int time_samples = 3000;
    float dt = 1e-3;
    float dx = 12.5;
    float dz = 12.5;
    float freq = 5.0;
    float a = M_PI*M_PI*freq*freq;

    float *wavelet = malloc(sizeof(float)*time_samples);;
    float *c = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *rho = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vx_0 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vx_1 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vx_2 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vz_0 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vz_1 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *vz_2 = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *Px = malloc(sizeof(float)*model_size[0]*model_size[1]);
    float *Pz = malloc(sizeof(float)*model_size[0]*model_size[1]);
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
        
        memcpy(temp_field, vx_0, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(vx_0, temp_field, model_size, dx, 'x');
        
        memcpy(temp_field, vz_0, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(vz_0, temp_field, model_size, dz, 'z');

        memcpy(P_2, P_1, sizeof(float)*model_size[0]*model_size[1]);
        memcpy(P_1, P_0, sizeof(float)*model_size[0]*model_size[1]);
        for(int x=0;x<model_size[0];x++){
            for(int z=0;z<model_size[1];z++){
                int idx = x+z*model_size[0];
                P_0[idx] = 2.0 * ( dt*rho[idx]*c[idx]*c[idx]*( vx_0[idx] + vz_0[idx] ) + P_2[idx]/2.0 );
            }
        }
        
        memcpy(vx_2, vx_1, sizeof(float)*model_size[0]*model_size[1]);
        memcpy(vx_1, vx_0, sizeof(float)*model_size[0]*model_size[1]);
        memcpy(temp_field, P_0, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(Px, temp_field, model_size, dx, 'x');
        for(int x=0;x<model_size[0];x++){
            for(int z=0;z<model_size[1];z++){
                int idx = x+z*model_size[0];
                vx_0[idx] = 2.0*dt*(-Px[idx] + vx_2[idx]/2.0)/rho[idx];
            }
        }

        memcpy(vz_2, vz_1, sizeof(float)*model_size[0]*model_size[1]);
        memcpy(vz_1, vz_0, sizeof(float)*model_size[0]*model_size[1]);
        memcpy(temp_field, P_0, sizeof(float)*model_size[0]*model_size[1]);
        spatial_diff(Pz, temp_field, model_size, dz, 'z');
        for(int x=0;x<model_size[0];x++){
            for(int z=0;z<model_size[1];z++){
                int idx = x+z*model_size[0];
                vz_0[idx] = 2.0*dt*(-Pz[idx] + vz_2[idx]/2.0)/rho[idx];
            }
        }

        P_0[source_loc[0] + model_size[0]*source_loc[1]] += wavelet[t];
        float max=0;
        for(int x=0;x<model_size[0]*model_size[1];x++){
            if (abs(P_0[x])>max){max=abs(P_0[x]);}
        }
        for(int x=0;x<model_size[0]*model_size[1];x++){
            temp_field[x] = 255*P_0[x]/max;
        }
        printf("max val is %f\n",max);
        // return -1;

        /* Lock the screen for direct access to the pixels */
        if ( SDL_MUSTLOCK(screenSurface) ) {
            if ( SDL_LockSurface(screenSurface) < 0 ) {
                fprintf(stderr, "Can't lock screen: %s\n", SDL_GetError());
                return -1;
            }
        }
        for(int x=0;x<model_size[0];x++){
            for(int z=0;z<model_size[1];z++){
                Uint32 Pcolor;
                int idx = x+z*model_size[0];
                Pcolor = SDL_MapRGB(screenSurface->format, (int)temp_field[idx], (int)temp_field[idx], (int)temp_field[idx]);
                putpixel(screenSurface, x, z, Pcolor);
            }
        }
        if ( SDL_MUSTLOCK(screenSurface) ) {
            SDL_UnlockSurface(screenSurface);
        }
        //Fill the surface white
        // SDL_FillRect( screenSurface, NULL, SDL_MapRGB( screenSurface->format, 0xFF, 0xFF, 0xFF ) );
        
        //Update the surface
        SDL_UpdateWindowSurface( window );

        //Hack to get window to stay up
        // SDL_Event e; bool quit = false; while( quit == false ){ while( SDL_PollEvent( &e ) ){ if( e.type == SDL_QUIT ) quit = true; } }

    }

    //Destroy window
    SDL_DestroyWindow( window );

    //Quit SDL subsystems
    SDL_Quit();

    free(c);
    free(rho);

    return 0;
}
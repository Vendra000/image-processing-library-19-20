/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

ip_mat * ip_mat_create (unsigned int h, unsigned int w, unsigned int k, float v){
    
    unsigned int i, j, l;
    
    ip_mat * result = (ip_mat *) malloc(sizeof(ip_mat));
    
    if(result != 0){
        result -> w = w; // larghezza, cioè numero di colonne
        result -> h = h; // altezza, cioè numero di righe
        result -> k = k; // profondità, cioè numero di canali
    
        result -> stat = (stats *)malloc(sizeof(stats)*k);
    
        result -> data = (float***) malloc(sizeof(float**) * h);
    
        for(i = 0; i < h; i++){
        
            result -> data[i] = (float**) malloc(sizeof(float*) * w);
        
            for(j = 0; j < w; j++){
            
                result-> data[i][j] = (float*) malloc(sizeof(float) * k);
            
                for(l = 0; l < k; l++){
                
                    set_val(result,i,j,l,v);
                
                }
            }
        }
    
    compute_stats(result);
    
    return result;
    
    }
    
    else {
        printf("errore create");
        exit(1);
    }
}

void ip_mat_free(ip_mat *a){
    
    if(a != 0){
        
        unsigned int i,j;
    
        free(a -> stat);
    
        for(i = 0;i < a->h; i++){
        
            for(j = 0;j < a->w; j++){
            
                free(a -> data[i][j]);
            }
        
            free(a -> data[i]);
        }
    
        free(a -> data);
        free(a);
    }
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

void compute_stats(ip_mat * t){
    
    float min, max, sum, mean;
    unsigned int i,j,l;
    
    for(l=0;l<t->k;l++){
        
        min = get_val(t,0,0,l);
        max = get_val(t,0,0,l);
        sum = 0.;
        
        for(i=0;i<t->h;i++){
            
            for(j=0;j<t->w;j++){
                
                if(get_val(t,i,j,l) < min)
                    min = get_val(t,i,j,l);
                
                if(get_val(t,i,j,l) > max)
                    max = get_val(t,i,j,l);
                
                sum += get_val(t,i,j,l);
            }
            
        }
        
        mean = sum / ((float)(t->h * t->w));
        
        t -> stat[l].min = min;
        t -> stat[l].max = max;
        t -> stat[l].mean = mean;
    }
}

void ip_mat_init_random(ip_mat *t, float mean, float var){
    
    unsigned int i,j,l;
    float rand;
    
    for(i=0;i<(t->h);i++){
        
        for(j=0;j<(t->w);j++){
            
            for(l=0;l<t->k;l++){
                
                rand = get_normal_random(mean, var);
                set_val(t,i,j,l,rand);
                
            }
        }
    }
    compute_stats(t);
}

ip_mat * ip_mat_copy(ip_mat * in){
    
    unsigned int i,j,l;
    ip_mat * new;
    
    new = ip_mat_create(in->h, in->w, in->k, 0.);
    
    for(i=0;i<(in->h);i++){
        
        for(j=0;j<(in->w);j++){
            
            for(l=0;l<(in->k);l++){
                
                set_val(new, i, j, l, get_val(in, i, j, l));
            }
        }
    }
    compute_stats(new);
    return new;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    
    unsigned int i,j,l;
    ip_mat * out = ip_mat_create(row_end-row_start+1, col_end-col_start+1, t->k, 0.);
    
    for(i=row_start;i<(row_end);i++){
        for(j=col_start;j<(col_end);j++){
            for(l=0;l<(t->k);l++){
                set_val(out,i-row_start,j-col_start,l,get_val(t,i,j,l));
            }
        }
    }
    
    return out;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    unsigned int i,j,l;
    ip_mat * out;
    
    if (dimensione == 0){
        
        out = ip_mat_create(a->h + b->h, a->w, a->k, 0.);
        
        for(i=0;i<(out->h);i++){
        
            for(j=0;j<(a->w);j++){
            
                for(l=0;l<(a->k);l++){
                
                    if (i<a->h)
                        set_val(out,i,j,l,get_val(a,i,j,l));
                    else
                        set_val(out,i,j,l,get_val(b,i-a->h,j,l));
                }
            }
        }
    }
    
    else if (dimensione == 1){
        
        out = ip_mat_create(a->h, a->w + b -> w, a->k, 0);
        
        for(i=0;i<(a->h);i++){
        
            for(j=0;j<(out->w);j++){
            
                for(l=0;l<(a->k);l++){
                
                    if (j<a->w)
                        set_val(out,i,j,l,get_val(a,i,j,l));
                    else
                        set_val(out,i,j,l,get_val(b,i,j-a->w,l));
                }
            }
        }
        
    }
    
    else if (dimensione == 2){
        
        out = ip_mat_create(a->h, a->w + b -> w, a->k, 0);
        
        for(i=0;i<(a->h);i++){
        
            for(j=0;j<(a->w);j++){
            
                for(l=0;l<(out->k);l++){
                
                    if (l<a->k)
                        set_val(out,i,j,l,get_val(a,i,j,l));
                    else
                        set_val(out,i,j,l,get_val(b,i,j,l-a->k));
                }
            }
        }
        
    }
    
    else{
        printf("Errore concat, dimensione non esistente !!!");
        exit(1);
    }
    
    return out;
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    
    ip_mat * out = ip_mat_create(a->h,a->w,a->k,0.);
    
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        
        unsigned int i,j,l;
        for(i=0;i<out->h;i++){
            for(j=0;j<out->w;j++){
                for(l=0;l<out->k;l++)
                    set_val(out,i,j,l,get_val(a,i,j,l)+get_val(b,i,j,l));
            }
        }
    }
    
    else{
        printf("dimensioni matrici sum non identiche");
        exit(1);
    }
    
    return out;
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    
    ip_mat * out = ip_mat_create(a->h,a->w,a->k,0.);
    
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        
        unsigned int i,j,l;
        
        for(i=0;i<out->h;i++){
            for(j=0;j<out->w;j++){
                for(l=0;l<out->k;l++)
                    set_val(out,i,j,l,get_val(a,i,j,l)-get_val(b,i,j,l));
            }
        }
    }
    
    else{
        printf("dimensioni matrici sub non identiche");
        exit(1);
    }
    
    return out;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    
    unsigned int i,j,l;
    ip_mat * out = ip_mat_copy(a);
    
    for(i=0;i<out->h;i++){
        for(j=0;j<out->w;j++){
            for(l=0;l<out->k;l++){
                set_val(out,i,j,l,get_val(out,i,j,l) * c);
            }
        }
    }
    
    return out;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    
    unsigned int i,j,l;
    ip_mat * out = ip_mat_copy(a);
    
    for(i=0;i<out->h;i++){
        for(j=0;j<out->w;j++){
            for(l=0;l<out->k;l++)
                set_val(out,i,j,l,get_val(out,i,j,l) + c);
        }
    }
    
    return out;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    
    ip_mat * out = ip_mat_copy(a);
    
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        
        unsigned int i,j,l;
        for(i=0;i<out->h;i++){
            for(j=0;j<out->w;j++){
                for(l=0;l<out->k;l++)
                    set_val(out,i,j,l,(get_val(a,i,j,l) + get_val(b,i,j,l))/2);
            }
        }
    }
    
    else{
        printf("dimensioni matrici non identiche");
        exit(1);
    }
    
    return out;
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    
    ip_mat * out = ip_mat_copy(in);
    unsigned int i,j,l;
    float sum,mean;
    
    for (i=0;i<out->h;i++){
        
        for(j=0;j<out->w;j++){
            
            sum = 0.;
            
            for(l=0;l<out->k;l++){
                sum += get_val(out,i,j,l);
                mean = sum /((float)(out -> k));
            }
            
            for (l=0;l<out->k;l++)
                set_val(out,i,j,l,mean);
        }
    }
    return out;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    
    ip_mat * img1 = ip_mat_copy(a);
    ip_mat * img2 = ip_mat_copy(b);
    
    ip_mat * out = ip_mat_sum( ip_mat_mul_scalar(img1, alpha), ip_mat_mul_scalar(img2,(1.0-alpha)) );
    
    ip_mat_free(img1);
    ip_mat_free(img2);
    
    return out;
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    ip_mat * out = ip_mat_copy(a);
    return ip_mat_add_scalar(out, bright);
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    ip_mat * crpt = ip_mat_copy(a);
    ip_mat * out;
    
    ip_mat_init_random(crpt,0.,amount);
    
    out = ip_mat_sum(a,crpt);
                
    return out;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0.);
    ip_mat * padded;
    unsigned int pad_h, pad_w, i, j, l, m, n;
    float acc = 0.;
    
    pad_h = ((f->h)-1)/2;
    pad_w = ((f->w)-1)/2;
    
    padded = ip_mat_padding(a, pad_h, pad_w);
    compute_stats(padded);
    
    for(l = 0; l < out -> k; l++){
        
        for(i = 0; i < out -> h; i++){
            
            for(j = 0; j < out -> w; j++){
                
                for(m = 0; m < f -> h; m++){
                    
                    for(n = 0; n < f -> w; n++){
                        
                        acc += get_val(f,m,n,l) * get_val(padded,m+i,n+j,l);
                    }
                }
                
                set_val(out,i,j,l,acc);
                acc = 0.;
                
            }
        }
    }
    
    ip_mat_free(padded);
    
    compute_stats(out);
    
    return out;
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
    unsigned int i,j,l;
    ip_mat * out = ip_mat_create(a->h + 2*pad_h, a->w + 2*pad_w, a->k, 0);
    
    for(i=0;i<a->h;i++){
        for(j=0;j<a->w;j++){
            for(l=0;l<a->k;l++){
                set_val(out,i + pad_h, j + pad_w, l, get_val(a,i,j,l));
            }
        }
    }
    
    compute_stats(out);
    
    return out;
}

ip_mat * create_sharpen_filter(){
    unsigned int i;
    ip_mat * out = ip_mat_create(3,3,3,0.);
    
    for(i=0;i<3;i++){
        set_val(out,0,1,i,-1.);
        set_val(out,1,0,i,-1.);
        set_val(out,1,2,i,-1.);
        set_val(out,2,1,i,-1.);
        set_val(out,1,1,i,5.);
    }
    return out;
}

ip_mat * create_edge_filter(){
    unsigned int i;
    ip_mat * out = ip_mat_create(3,3,3,-1.);
    
    for(i=0;i<3;i++){
        set_val(out,1,1,i,8.);
    }
    return out; 
}

ip_mat * create_emboss_filter(){
    unsigned int i;
    ip_mat * out = ip_mat_create(3,3,3,1.);
    
    for(i=0;i<3;i++){
        set_val(out,0,0,i,-2.);
        set_val(out,0,1,i,-1.);
        set_val(out,0,2,i,0.);
        set_val(out,1,0,i,-1.);
        set_val(out,2,0,i,0.);
        set_val(out,2,2,i,2.);
    }
    return out;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
    ip_mat * out;
    float c;
    
    c = 1.0/(w*h);
    out = ip_mat_create(h,w,k,c);
    
    return out;
}

ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma){
    
    int cx, cy, x, y;
    unsigned int i, j, l;
    float g, sigma2, sum;
    ip_mat * out = ip_mat_create(h,w,k,0.);
    cx = (w-1)/2;
    cy = (h-1)/2;
    sigma2 = pow(sigma,2);
    
    for(l=0;l<k;l++){
        
        sum = 0.;
        
        for(i=0;i<h;i++){
            for(j=0;j<w;j++){
                x = i-cx;
                y = j-cy;
                g = (1./(2.*PI*sigma2)) * exp( (-((pow(x,2)) + (pow(y,2)))/(2*sigma2) ));
                sum += g;
                set_val(out,i,j,l,g);
            }
        }
        
        for(i=0;i<h;i++){
            for(j=0;j<w;j++){
                set_val(out,i,j,l,get_val(out,i,j,l)/sum);
            }
        }
    }
    
    return out;
}

void rescale(ip_mat * t, float new_max){
    
    unsigned int i,j,l;
    float min, max, new_val;
    compute_stats(t);
    
    for(l=0; l<t->k; l++){
        min = t -> stat[l].min;
        max = t -> stat[l].max;
        for(i=0; i<t->h; i++){
            for(j=0; j<t->w; j++){
                new_val = ((get_val(t,i,j,l)) - min) / (max-min);
                set_val(t,i,j,l,new_val * new_max);
            }
        }
    } 
}

void clamp(ip_mat * t, float low, float high){
    
    unsigned int i,j,l;
    
    for(i=0;i<t->h;i++){
        for(j=0;j<t->w;j++){
            for(l=0;l<t->k;l++){
                
                if(get_val(t,i,j,l) < low){
                    set_val(t,i,j,l,low);
                }
                
                if(get_val(t,i,j,l) > high){
                    set_val(t,i,j,l,high);
                }
            }
        }
    }
}



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define H_STR "H"
#define F_PRIME "f_prime"
#define INITIAL_STR "initial_coord"
#define DIST_STR "distance"
#define PERTURBE_STR "BEFORE_CG"
#define FINAL_STR "AFTER_CG"
#define NEIGH_STR "neigh"

#define UX 8 //(4 x 2)
#define UY 8
#define UZ 8
#define A_CONSTANT 5.64/2
#define SMALL_SIGMA 1E-3

#define A_Na_Na 7895.4
#define RHO_Na_Na 0.1709
#define C_Na_Na 29.06

#define A_Na_Cl 2314.7
#define RHO_Na_Cl 0.2903
#define C_Na_Cl 0.

#define A_Cl_Cl 1227.2
#define RHO_Cl_Cl 0.3214
#define C_Cl_Cl  29.06

#define Z_NA 0.988
#define Z_CL 0.988

#define COULOMB (Z_NA * Z_CL)/(4.*M_PI*0.005526349406) 

double distance_cutoff = 10 * 5.64;  
int atom_num = UX * UY * UZ;
double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors
int i_max = 2;

void FILE_WRITING(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LLAB4_2_%s.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \t %lf \t %lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2], arr[i][3]);
    }
    printf("%s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

void DIST_FILE_WRITING(char *lattice_structure, double** position_vector , double* distance_arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LAB4_%s_neighb_list.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
       fprintf(filepointer,"%lf \t %lf \t %lf \t %lf\n",position_vector[i][0], position_vector[i][1] , position_vector[i][2], distance_arr[i]);
    }
    printf("distance list of %s has been printed!\n",lattice_structure);
    fclose(filepointer);
}


void DEBUG_FILE_WRITING(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LAB4_DEBUG_%s.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \t %lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }
    printf("%s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

void CG_FILE_WRITING(double *arr, int num){    
    FILE *filepointer = fopen("LAB4_2_CG_energy_vs_steps.txt", "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%0.8lf \n",arr[i]);
    }
    printf("energy vs steps has been printed!\n");
    fclose(filepointer);
}

void neigh_file_writing(char *lattice_structure, double *neigh_dist, int** neigh_list, double* interaction_list, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LAB4_ion_%s_list.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
       fprintf(filepointer,"%d \t %d \t %lf \t %lf \n",neigh_list[i][0], neigh_list[i][1] , neigh_dist[i], interaction_list[i]);
    }
    printf("neighbour list of %s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

/*b1 b2 b3 VECTOR GENERATION=============================================================================================*/
void fn_a_vect(){
    //generating a1 a2 a3 vectors
    a[0][0] = UX * A_CONSTANT;
    a[1][1] = UY * A_CONSTANT;
    a[2][2] = UZ * A_CONSTANT;
    for(int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < 3 ; j++){
            if(i != j){
                a[i][j] = 0;
            }
        }
    }
    return;
}

void fn_cross_product(int idx1 , int idx2 , int idx3){ //b_(idx1 + 1) = a_(idx2 + 1) x a_(idx3 + 1)
        b[idx1][0] = a[idx2][1]*a[idx3][2] - a[idx2][2]*a[idx3][1];
        b[idx1][1] = a[idx2][2]*a[idx3][0] - a[idx2][0]*a[idx3][2];
        b[idx1][2] = a[idx2][0]*a[idx3][1] - a[idx2][1]*a[idx3][0];
    return;
}

double fn_2pi_triple_product(){
    volume = (a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2]);
    double answer = 1/volume;
    return answer;
}

void fn_reciprocal(){
    fn_cross_product(0 , 1 , 2); // b1 = a2 x a3
    double scalars = fn_2pi_triple_product(); //a1 dot (a2 x a3)
    fn_cross_product(1 , 2 , 0); // b2 = a3 x a1
    fn_cross_product(2 , 0 , 1); // b3 = a1 x a2

    for (int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < 3 ; j++){
            b[i][j] = b[i][j] * scalars;
        }
    }
    return;
}

/*ATOM COORDS GENERATION==============================================================================================*/

void fn_simplecubic(double **coords){
    int atom_type = 1;
    for (int i = 0; i < atom_num ; i++) {
        coords[i][0] = (i % UX) * A_CONSTANT;
        coords[i][1] = (i / UX) % UY * A_CONSTANT;
        coords[i][2] = (i / UX / UY) % UZ * A_CONSTANT;

        atom_type *= -1;
        if(UX % 2 == 0 && i % (UX * UY) == 0){
            atom_type *= -1;
        }
        if(UX % 2 == 0 && i % UX == 0){ //if perodicity of y is an even number and the loop has reached to the end, then the atom in the next y 
            atom_type *= -1;
        }
        coords[i][3] = atom_type; //all the 'simple cubic' lattice points are Na
  
    }
    return;
}

/*NEIGHBOUR LIST GENERATION=======================================================================================*/
double fn_pbc(double* t){
    double n[3];
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2]) , 1); //generate fractional numbers
        //apply periodic boundary conditions [-0.5 , 0.5)
        if(n[i] > 0.5){
            difference_plus = n[i] - 0.5;
            n[i] = - 0.5 + difference_plus;
            difference_plus = 0;
        }
        if(n[i] <= -0.5){
            difference_minus = -n[i] - 0.5;
            n[i] = 0.5 - difference_minus;
            difference_minus = 0;
        }
    }
    for(int i = 0 ; i < 3 ; i++){
        t[i] = n[0] * a[0][i] + n[1] * a[1][i] + n[2] * a[2][i]; //fractional coordinate
    }
    double distance = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    return distance;
}


int fn_neighbour_list(int **neigh_list, double **atom_coords, int num, double **vectors,double* neigh_distance, double *interaction_list){
    double distance = 0;
    int nearest_neighbour_num = 0;
    double t[3];
    double interaction_type = 0; //if both are Na atoms, interaction = -2, if one is Na , one is Cl, interaction = 0 , if both are Cl, interaction = 2


    for(int i = 0 ; i < num ; i++){
        for(int j = 0 ; j < num ; j++){ //j = i + 1 to avoid double counting (1 2 , 2 1), i + 1 for skipping (0 0) , (1 1)
            if(i == j){continue;}

            for(int k = 0 ; k < 3 ; k++){
                t[k] = atom_coords[j][k] - atom_coords[i][k];
            }
            interaction_type = atom_coords[j][3] + atom_coords[i][3];

            distance = fn_pbc(t);

            if(distance <= distance_cutoff){ //evaulate whether atom i and atom j are nearest neighbours
                for(int k = 0; k < 3 ; k++){
                    vectors[nearest_neighbour_num][k] = t[k];
                }

                neigh_list[nearest_neighbour_num][0] = i;
                neigh_list[nearest_neighbour_num][1] = j;
                
                interaction_list[nearest_neighbour_num] = interaction_type;

                neigh_distance[nearest_neighbour_num] = distance;

                nearest_neighbour_num = nearest_neighbour_num + 1;
            }
        }
    }

    return nearest_neighbour_num;
}


double fn_coulomb_buck_potential(int neigh_size, double* neigh_distance , double* interaction_list){
    double potential = 0;
    double temp_potential = 0;
    double r_6 = 0;
    double r = 0; 

    for(int i = 0 ; i < neigh_size ; i++){
        r = neigh_distance[i];
        r_6 = r * r * r * r * r * r;
        if(interaction_list[i] == -2){ //Na-Na interaction
            temp_potential = COULOMB / r + A_Na_Na * exp(-r /RHO_Na_Na) - C_Na_Na / r_6 ;
        }
        if(interaction_list[i] == 0){ //Na-Cl interaction
            temp_potential = (-COULOMB  / r)+ (A_Na_Cl * exp(-r / RHO_Na_Cl)) - (C_Na_Cl / r_6);
        }
        if(interaction_list[i] == 2){ //Cl-Cl interaction
            temp_potential = ( COULOMB / r) + A_Cl_Cl * exp(-r / RHO_Cl_Cl) - (C_Cl_Cl / r_6);
        }
        potential += temp_potential;
    }

    //return 2* potential;
    return potential/atom_num;
}

/*INITIALISATION=====================================================================================================*/
void FN_initialise(double** atom_coords){

    fn_a_vect(); //initialise a1 a2 a3 cell vectors
    fn_reciprocal(); //calculate the reciprocal vectors and the volume
    fn_simplecubic(atom_coords); //generate coordinates

    return;
}


/*PERTURBATING ATOMS================================================================================================*/
void FN_perturbate(double** atom_coords){
    double min = -0.025;
    double max = 0.3;
    for(int i = 0 ; i < atom_num ; i++){
        for(int j = 0 ; j < 3 ; j++){
            double random = (rand() * (0.1) / RAND_MAX ) - 0.05;
            atom_coords[i][j] += A_CONSTANT * random;
        }
    }

    return;
}


double fn_F_prime_coulomb(double r,double interaction_type, double vector_component){
    double F_prime = 0;
    double r3 = r * r * r;


        if(interaction_type == -2){ //Na-Na interaction
            F_prime = (- COULOMB / r3);
        }
        if(interaction_type == 0){ //Na-Cl interaction
            F_prime = ( - COULOMB / r3);
        }
        if(interaction_type == 2){ //Cl-Cl interaction
            F_prime = (- COULOMB / r3);
        }

    return  F_prime * vector_component;
}

double fn_F_prime_buck(double r,double interaction_type,double vector_component){
    double F_prime = 0;
    double r8 = r * r * r * r * r * r * r * r ;

        if(interaction_type == -2){ //Na-Na interaction
            F_prime =  - (A_Na_Na / RHO_Na_Na * exp(-r / RHO_Na_Na)) + (6 * C_Na_Na / r8);
        }
        if(interaction_type == 0){ //Na-Cl interaction
            F_prime =  - ( A_Na_Cl / RHO_Na_Cl * exp(-r / RHO_Na_Cl)) + (6 * C_Na_Cl / r8);
        }
        if(interaction_type == 2){ //Cl-Cl interaction
            F_prime =  - ( A_Cl_Cl / RHO_Cl_Cl * exp(-r / RHO_Cl_Cl)) + (6 * C_Cl_Cl / r8);
        }

    return  F_prime*vector_component;
}

void fn_coulomb_buck_derivative(int** neigh_list,int neigh_size, double* neigh_distance , double* interaction_list, double** vector, double** gradient, int idx0 , int idx1){

    for(int i = 0 ; i < neigh_size ; i++){
        //calculate the magnitude of F prime
        for(int j = 0 ; j < 3 ; j++){
            double F_prime =  (fn_F_prime_coulomb(neigh_distance[i],interaction_list[i],vector[i][j]) + fn_F_prime_buck(neigh_distance[i],interaction_list[i],vector[i][j]));
            //gradient of atom A in the pair consists of x y z component
            gradient[neigh_list[i][idx0]][j] += -  F_prime/2 ;
            //gradient of atom B in the pair also have x y z component, but its gradient is the negative of the gradient of atom A
            gradient[neigh_list[i][idx1]][j] +=  - F_prime/2;
        }

    }



    return;
}

void fn_coulomb_buck_derivative2(int** neigh_list,int neigh_size, double* neigh_distance , double* interaction_list, double** vector, double** gradient, int idx0 , int idx1){

    for(int i = 0 ; i < neigh_size ; i++){
        //calculate the magnitude of F prime


        for(int j = 0 ; j < 3 ; j++){
            double F_prime =  (fn_F_prime_coulomb(neigh_distance[i],interaction_list[i],vector[i][j]) + fn_F_prime_buck(neigh_distance[i],interaction_list[i],vector[i][j]));

            //gradient of atom A in the pair consists of x y z component
            gradient[neigh_list[i][idx0]][j] +=   F_prime /2;
            //gradient of atom B in the pair also have x y z component, but its gradient is the negative of the gradient of atom A
            gradient[neigh_list[i][idx1]][j] +=   F_prime/2;
        }
        /*double F_prime =  fn_F_prime_coulomb(neigh_distance[i],interaction_list[i]) + fn_F_prime_buck(neigh_distance[i],interaction_list[i]);

        for(int j = 0 ; j < 3 ; j++){
            //gradient of atom A in the pair consists of x y z component
            gradient[neigh_list[i][idx0]][j] +=  F_prime * vector[i][j];
        }

        for(int j = 0 ; j < 3 ; j++){
            //gradient of atom B in the pair also have x y z component, but its gradient is the negative of the gradient of atom A
            gradient[neigh_list[i][idx1]][j] -=  F_prime * vector[i][j];
        }*/ 

    }



    return;
}

double fn_line_minimisation(double** atom_coords,double** H,double** f_prime,int step){
     /*INITIALISE ARRAYS FOR ATOMS WHICH WILL BE SLIGHTLY DIPLACED FOR LINE MINIMISATION===============================*/

    double** little_displacement_atom_coords = (double **)malloc(atom_num * sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        little_displacement_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    double row_num = atom_num * (atom_num - 1);

    int **neigh_list = (int **)malloc(row_num * sizeof(int *));
    for (int i = 0; i < row_num; i++){
        neigh_list[i] = (int *)malloc(2 * sizeof(int));
    }

    double **vector = (double **)malloc(row_num * sizeof(double *));
    for (int i = 0; i < row_num; i++){
        vector[i] = (double *)malloc(3 * sizeof(double));
    }

    double *neigh_distance = (double*)malloc(row_num * sizeof(double));

    double **displaced_f_prime = (double **)malloc(atom_num * sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        displaced_f_prime[i] = (double *)malloc(3 * sizeof(double));
    }

    double *displaced_interaction_list = (double*)malloc(row_num * sizeof(double));

    //atom coords are shifted slightly
    for(int i = 0 ; i < atom_num ; i++){
        for(int j = 0 ; j < 3 ; j++){
            little_displacement_atom_coords[i][j]  = (SMALL_SIGMA * H[i][j]) + atom_coords[i][j];
        }
    }

    //create new neighbour list for the displaced atoms
    int displaced_neigh_num = fn_neighbour_list(neigh_list,little_displacement_atom_coords, atom_num, vector, neigh_distance,displaced_interaction_list);

    //calculate displaced f prime
    fn_coulomb_buck_derivative2(neigh_list,displaced_neigh_num,neigh_distance,displaced_interaction_list,vector,displaced_f_prime,0,1);

    //find numerator & denominator
    double numerator = 0 ; 
    double denominator = 0;
    double alpha = 0;
    for(int i = 0 ; i < atom_num ; i++){
        numerator += (f_prime[i][0] * H[i][0]) + (f_prime[i][1] * H[i][1]) + (f_prime[i][2] * H[i][2]);
        denominator += ((displaced_f_prime[i][0] - f_prime[i][0]) * H[i][0]) + ((displaced_f_prime[i][1] - f_prime[i][1]) * H[i][1]) + ((displaced_f_prime[i][2] - f_prime[i][2]) * H[i][2]);
    }

    if(step == 0){
        DEBUG_FILE_WRITING(H_STR,H,atom_num);
        DEBUG_FILE_WRITING(F_PRIME,f_prime,atom_num);
    }

    alpha =  -SMALL_SIGMA * numerator / (denominator);

    printf("numerator = %lf \t denominator = %lf \t alpha = %lf \n", numerator, denominator, alpha);



    free(little_displacement_atom_coords);
    free(displaced_f_prime);
    free(neigh_distance);
    free(vector);
    free(neigh_list);
    free(displaced_interaction_list);

    return alpha;
}

double fn_gamma(double** g1 , double** g2){
    double numerator = 0;
    double denominator = 0;
    for(int i = 0 ; i < atom_num ; i++){
        numerator += (g2[i][0] * g2[i][0]) + (g2[i][1] * g2[i][1]) + (g2[i][2] * g2[i][2]);
        denominator += (g1[i][0] * g1[i][0]) + (g1[i][1] * g1[i][1]) + (g1[i][2] * g1[i][2]);
    }

    return numerator / denominator;
}

int FN_CG(double** atom_coords,double* coulomb_energy){
    double residual_magnitude = 1;
    int i = 0;

    /*INITIALISE ARRAYS============================================================================
        
        neighbour_list = stores neighbour pair atom labels

        position_vectors = stores the distance vector of the atom pair

        neighbour_distance = the distance between the two atom pair

        interaction_list = stores whether the atom pair is Na-Na, Na-Cl, or Cl-Cl
        
    */
    
    double row_num = atom_num * (atom_num - 1);

    int **neighbour_list = (int **)malloc(row_num * sizeof(int *));
        for (int j = 0; j < row_num; j++){
            neighbour_list[j] = (int *)malloc(2 * sizeof(int));
        }

    double **position_vector = (double **)malloc(row_num * sizeof(double *));
        for (int j = 0; j < row_num; j++){
            position_vector[j] = (double *)malloc(3 * sizeof(double));
        }

    double *neighbour_distance = (double*)malloc(row_num * sizeof(double));

    double *interaction_list = (double*)malloc(row_num * sizeof(double));

    double **g1 = (double **)malloc(atom_num * sizeof(double *));
        for (int j = 0; j < atom_num; j++){
            g1[j] = (double *)malloc(3 * sizeof(double));
        }

    double **h = (double **)malloc(atom_num * sizeof(double *));
        for (int j = 0; j < atom_num; j++){
            h[j] = (double *)malloc(3 * sizeof(double));
        }


    //obtain the neighbour list
    int neighbour_num = 0;
    neighbour_num = fn_neighbour_list(neighbour_list, atom_coords, atom_num , position_vector, neighbour_distance,interaction_list);
    neigh_file_writing(NEIGH_STR , neighbour_distance , neighbour_list,interaction_list , neighbour_num);
    DIST_FILE_WRITING(DIST_STR,position_vector,neighbour_distance,neighbour_num);

    //obtain g1 by calculating the coulomb buck derivative at all atom coordinates
    fn_coulomb_buck_derivative(neighbour_list,neighbour_num,neighbour_distance, interaction_list, position_vector,g1,0,1);
    //obtain h by calculating the coulomb buck derivative at all atom coordinates
    fn_coulomb_buck_derivative(neighbour_list,neighbour_num,neighbour_distance,interaction_list, position_vector,h,0,1);
        

    while(i < i_max || residual_magnitude > 0.05){

        int **this_neighbour_list = (int **)malloc(row_num * sizeof(int *));
            for (int j = 0; j < row_num; j++){
                this_neighbour_list[j] = (int *)malloc(2 * sizeof(int));
            }
        double **f_prime = (double **)malloc(atom_num * sizeof(double *));
        for (int j = 0; j < atom_num; j++){
            f_prime[j] = (double *)malloc(3 * sizeof(double));
        }
        double** this_atom_coords = (double **)malloc(atom_num * sizeof(double *));
        for (int i = 0; i < atom_num; i++){
            this_atom_coords[i] = (double *)malloc(4 * sizeof(double));
        }
        double **g2 = (double **)malloc(atom_num * sizeof(double *));
        for (int j = 0; j < atom_num; j++){
            g2[j] = (double *)malloc(3 * sizeof(double));
        }
        double *this_neighbour_distance = (double*)malloc(row_num * sizeof(double));

        double *this_interaction_list = (double*)malloc(row_num * sizeof(double));

        double **this_position_vector = (double **)malloc(row_num * sizeof(double *));
        for (int j = 0; j < row_num; j++){
            this_position_vector[j] = (double *)malloc(3 * sizeof(double));
        }



        coulomb_energy[i] = fn_coulomb_buck_potential(neighbour_num,neighbour_distance, interaction_list);
            //store the coulomb energy at this CG step

        printf("coulomb energy = %lf \n", coulomb_energy[i]);

        int f_prime_neighbour_num = fn_neighbour_list(neighbour_list, atom_coords, atom_num , position_vector, neighbour_distance,interaction_list);
        //calculate f prime
        fn_coulomb_buck_derivative2(neighbour_list,f_prime_neighbour_num,neighbour_distance,interaction_list, position_vector,f_prime,0,1);

        /*PERFORM A SINGLE CG STEP====================================================================*/



        //obtain gamma through line minimisation
        double alpha = fn_line_minimisation(atom_coords,h,f_prime,i);

        //shift all atom coordinates
        for(int j = 0; j < atom_num ; j++){
            for(int k = 0 ; k < 3 ; k++){
                this_atom_coords[j][k] = atom_coords[j][k] + alpha * h[j][k];
            }
        }

        /*CALCULATE RESIDUAL MAGNITUDE TO CHECK IF ANOTHER CG STEP SHOULD BE TAKEN=====================*/
        residual_magnitude = 0;
        for(int j = 0; j < atom_num ; j++){
            residual_magnitude += (alpha * alpha) * (g1[j][0]*g1[j][0] +  g1[j][1]*g1[j][1] +  g1[j][2]*g1[j][2]);
    
        }

        //find neighbour number of the shifted atoms
        int this_neighbour_num = fn_neighbour_list(this_neighbour_list, this_atom_coords, atom_num , this_position_vector, this_neighbour_distance,this_interaction_list);

        //obtain g2 by calculating the coulomb buck derivative at all atom coordinates
        fn_coulomb_buck_derivative(this_neighbour_list,this_neighbour_num,this_neighbour_distance,this_interaction_list, this_position_vector,g2,0,1);



        double gamma = fn_gamma(g1,g2);

        //DEBUG_FILE_WRITING(H_STR,g2,atom_num);

        printf("gamma = %lf \t resi mag = %e \t alpha = %lf\n",gamma, residual_magnitude,alpha);

        for(int j = 0; j < atom_num ; j++){
            for(int k = 0 ; k < 3 ; k++){
                h[j][k] = g2[j][k] + (gamma * h[j][k]);
                g1[j][k] = g2[j][k];
                atom_coords[j][k] = this_atom_coords[j][k];
            }
        }

        for(int j = 0; j < this_neighbour_num ; j++){
            neighbour_distance[j] = this_neighbour_distance[j];
            interaction_list[j] = this_interaction_list[j];
            for(int k = 0 ; k < 2 ; k++){
                neighbour_list[j][k] = this_neighbour_list[j][k];
            }
            for(int k = 0 ; k < 3 ; k++){
                position_vector[j][k] = this_position_vector[j][k];
            }
        }

        neighbour_num = this_neighbour_num;

        free(f_prime);
        free(this_atom_coords);
        free(g2);
        free(this_interaction_list);
        free(this_neighbour_distance);
        free(this_neighbour_list);
        free(this_position_vector);
        i++;

    }

    FILE_WRITING(FINAL_STR ,atom_coords,atom_num);

    free(neighbour_list);
    free(neighbour_distance);
    free(position_vector);
    free(interaction_list);
    free(g1);
    free(h);

    return i;
}

/*MAIN==============================================================================================================*/
int main(){
    srand(time(NULL));

    double** fcc_atom_coords = (double **)malloc(atom_num * sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        fcc_atom_coords[i] = (double *)malloc(4 * sizeof(double));
    }
    //create array to hold energy for each SD step
    double *coulomb_energy;
    coulomb_energy = (double*)malloc(i_max * sizeof(double));

    FN_initialise(fcc_atom_coords);
    FILE_WRITING(INITIAL_STR , fcc_atom_coords,atom_num);
    //initial unpreturbed atom coords

    //perturbed atom coords
    FN_perturbate(fcc_atom_coords);
    //FILE_WRITING(PERTURBE_STR , fcc_atom_coords,atom_num);

    //perform steepest descent
    int CG_steps = FN_CG(fcc_atom_coords,coulomb_energy);
    CG_FILE_WRITING(coulomb_energy,CG_steps);

    free(fcc_atom_coords);
    free(coulomb_energy);
    
    return 0;
}
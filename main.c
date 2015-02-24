#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define NUM_COLORS 256

typedef struct Complex {
	double real;
	double imag;
}Complex;

// functie ce aduna 2 numere complexe
Complex add(double real1, double imag1, double real2, double imag2){
	Complex res;
	res.real = real1 + real2;
	res.imag = imag1 + imag2;
	return res;
}

// functie ce inmulteste 2 numere complexe
Complex multiply(double real1, double imag1, double real2, double imag2){
	Complex res;
	res.real = real1 * real2 - imag1 * imag2;
	res.imag = real1 * imag2 + imag1 * real2;
	return res;
}

// functie ce calculeaza modulul unui numar complex
double abs_value(double real, double imag){
	return sqrt(real * real + imag * imag);
}

int main(int argc, char *argv[]) {
	int  rank, type, max_steps, size;
	double x_min, x_max, y_min, y_max, resolution, x, y;	
        int width, height, i, j;
	FILE *fr, *fw;
	MPI_Init(&argc,&argv);
	MPI_Status status; 
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	// daca este masterul, citeste din fisier datele de intrare
	if(rank == 0){ 
		fr = fopen(argv[1],"r");
		fscanf(fr,"%d", &type);
		if(type == 0)  		
			fscanf(fr,"%lf %lf %lf %lf %lf %d",
				&x_min, &x_max, &y_min, &y_max, &resolution, &max_steps);
		else
			fscanf(fr,"%lf %lf %lf %lf %lf %d %lf %lf",
				&x_min, &x_max, &y_min, &y_max, &resolution, &max_steps, &x, &y);
		fclose(fr);
		// apoi trimite fiecarui proces datele din fisierul de intrare
		for (i = 1; i < size; i++) {
			MPI_Send(&type, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&x_min, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&x_max, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&y_min, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&y_max, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&resolution, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			MPI_Send(&max_steps, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			if (type == 1) {
				MPI_Send(&x, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				MPI_Send(&y, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			}
		}
	}
	else {
		//celelalte procese primesc datele de intrare de la master
		MPI_Recv(&type, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&x_min, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&x_max, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_min, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_max, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&resolution, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&max_steps, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		if (type == 1) {
			MPI_Recv(&x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&y, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
		}	
 	}
 	// se calculeaza latimea si inaltimea dupa formule
 	width = (int) floor((x_max-x_min)/resolution);
 	height = (int) floor((y_max-y_min)/resolution);
   	double z_real, z_imag;
   	int step, pixel, k;
   	int data = 1;
   	int* pixels = (int*) malloc ((width/size+ 1) * sizeof(int));
   	// se implemeanteaza algoritmul Julia
   	if(type == 1){
   		double z_copy_real, z_copy_imag; 
   		for(z_imag = y_min; z_imag < y_max; z_imag += resolution){
   			// voi trimite masterului cate o linie cu width/size componente
   			// pe prima pozitie voi retine linia pe care o trimit
   			int line = (int)floor((z_imag - y_min) / resolution);
			k = 0;
			pixels[k] = line;
			k++;
			// impart latimea pentru ca fiecare proces sa execute o parte
   			for( z_real =  x_min + rank * ((x_max-x_min)/size);
   			 z_real < x_min + (rank + 1) * ((x_max-x_min)/size); z_real += resolution){
				step = 0;
				z_copy_real = z_real;
				z_copy_imag = z_imag;
				while(abs_value(z_copy_real, z_copy_imag) < 2.0 && step < max_steps){
					double prod_real, prod_imag;
					prod_real = multiply(z_copy_real, z_copy_imag, z_copy_real, z_copy_imag).real;
					prod_imag = multiply(z_copy_real, z_copy_imag, z_copy_real, z_copy_imag).imag;
					z_copy_real = add(prod_real, prod_imag, x, y).real;
					z_copy_imag = add(prod_real, prod_imag, x, y).imag;
					step = step + 1;
				}
				pixel = step % NUM_COLORS;
				// pun pixelul calculat in vector
				pixels[k] = pixel;	
				k++;
			}
			// trimit un mesaj prin care masterul
			// sa stie ca o sa trimit vectori cu elementele calculate
			MPI_Send(&data, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			// trimit vectorul cu valorile calculate
			MPI_Send(pixels, width/size + 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}
   	}
   	
   	// se implemeanteaza algoritmul Mandelbrot
   	else if(type == 0){
   		double c_real, c_imag;
   		for(c_imag = y_min; c_imag < y_max; c_imag += resolution){
   			// voi trimite masterului cate o linie cu width/size componente
   			// pe prima pozitie voi retine linia pe care o trimit
   			int line = (int) floor((c_imag - y_min) / resolution);
			k = 0;
			pixels[k] = line;
			k++;
   			for( c_real =  x_min + rank * ((x_max-x_min)/size);
   			 c_real < x_min + (rank + 1) * ((x_max-x_min)/size); c_real += resolution){
				step = 0;
				z_real = 0;
				z_imag = 0;
				while(abs_value(z_real, z_imag) < 2.0 && step < max_steps){
					double prod_real, prod_imag;
					prod_real = multiply(z_real, z_imag, z_real, z_imag).real;
					prod_imag = multiply(z_real, z_imag, z_real, z_imag).imag;
					z_real = add(prod_real, prod_imag, c_real, c_imag).real;
					z_imag = add(prod_real, prod_imag, c_real, c_imag).imag;
					step = step + 1;
				}
				pixel = step % NUM_COLORS;
				// pun pixelul calculat in vector
				pixels[k] = pixel;	
				k++;
			}
			// trimit un mesaj prin care masterul sa stie 
			// ca o sa trimit vectori cu elementele calculate
			MPI_Send(&data, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			// trimit vectorul cu valorile calculate
			MPI_Send(pixels, width/size + 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}
   	}
   	// trimit un mesaj prin care masterul sa stie 
   	// ca am terminat de trimis valorile calculate 
   	data = 0;
	MPI_Send(&data, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	// masterul va scrie in fisierul de iesire datele calculate
   	if(rank == 0){
   		fw = fopen(argv[2],"w");
   		fprintf(fw, "P2\n");
   		fprintf(fw, "%d %d\n", width, height);
   		fprintf(fw, "%d\n", NUM_COLORS - 1);
   		// aloc memorie pentru matrice
   		int **mat = (int **)malloc(height * sizeof(int*));
   		for (i = 0; i < height; i++) {
			mat[i] = (int*) malloc (width * sizeof(int));
		}
		// masterul primeste de la fiecare proces in parte
		// vectori cu datele calculate
   		for (i = 0; i < size; i++) {
			MPI_Recv(&data, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			// cat timp se mai primesc valori
			while (data == 1) {
				MPI_Recv(pixels, width/size + 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
				// pe prima pozitie se afla pozitia liniei
				int line = pixels[0];
				// coloanele se afla la pozitia (k - 1) + i * width/size,
				// in functie de procesul care se executa
				for (k = 1 ; k <= width / size; k++) {
					mat[line][(k - 1) + i * width/size] = pixels[k];
				}
			// primesc mesaj sa stiu daca mai primesc valori sau nu
			MPI_Recv(&data, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			}
			
   		}
   		// scriu matricea in fisier, tinand cont 
   		// sa scriem coordonatele Y in ordine inversa
   		for (i = height - 1; i > 0; i--) {
			for (j = 0; j < width; j++) {
				fprintf(fw, "%d ", mat[i][j]);
			}
			fprintf(fw, "\n");
		}
   		fclose(fw);
   		// dezaloc memoria folosita
   		free(pixels);
   		for (i = 0; i < height; i++)
   			free(mat[i]);
   		free(mat);
   	}
      	MPI_Finalize();
}

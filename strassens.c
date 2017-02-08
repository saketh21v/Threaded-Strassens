#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void prnMtrx(int** a, int n, int m){
    printf("\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            printf("%d\t", a[i][j]);
        }
        printf("\n");
    }
}

int** createMatrix(int n){
    int** m = (int**) malloc(sizeof(int*) * n);
    for(int i=0;i<n;i++){
        m[i] = (int*) malloc(sizeof(int) * n);
    }
    return m;
}

int** createMtrx(int n, int m){
    int** a = (int**) malloc(sizeof(int*) * n);
    for(int i=0;i<n;i++){
        a[i] = (int*) malloc(sizeof(int) * m);
    }
    return a;
}


void deleteMtrx(int** a, int n){
    for(int i=0;i<n;i++){
        free(a[i]);
    }
    free(a);
}

void mtrxMul(int** a, int** b, int** r, int n){
    int sum = 0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                sum = sum + a[i][k]*b[k][j];
            }
            r[i][j] = sum;
            sum = 0;
        }
    }
}

void addMatrices(int** m1, int** m2, int** result, int a1, int a2, int b1, int b2, int n, int sub){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            if(sub == 0)
                result[i][j] = m1[a1 + i][a2 + j] + m2[b1 + i][b2 + j];
            else
                result[i][j] = m1[a1 + i][a2 + j] - m2[b1 + i][b2 + j];
    }
}

void addMtrx(int** a, int** b, int** r, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            r[i][j] = a[i][j] + b[i][j];
        }
    }
}


void mtrxMulRecur(int** a, int** b, int** r, int a1, int a2, int b1, int b2,int n){
    /*
        If n > 2
        Divide a into a11, a12, a21, a22 and b into b11, b12, b21, b22
        r11 = a11 * b11 + a12 * b21
        r12 = a11 * b12 + a12 * b22
        r21 = a21 * b11 + a22 * b21
        r22 = a21 * b12 + a22 * b22

        Let
            a11b11 = a11 * b11
            a12b21 = a12 * b21
            a11b12 = a11 * b12
            a12b22 = a12 * b22
            a21b11 = a21 * b11
            a22b21 = a22 * b21
            a21b12 = a21 * b12chr
            a22b22 = a22 * b22
    */
    int sum = 0;
    if(n > 2){
        int n2 = n/2;

        int** a11b11 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a11b11,a1,a2,b1,b2,n2);

        int** a12b21 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a12b21,a1,a2+n2,b1+n2,b2,n2);

        int** a11b12 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a11b12,a1,a2,b1,b2+n2,n2);

        int** a12b22 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a12b22,a1,a2+n2,b1+n2,b2+n2,n2);

        int** a21b11 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a21b11,a1+n2,a2+0,b1+0,b2+0,n2);

        int** a22b21 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a22b21,a1+n2,a2+n2,b1+n2,b2+0,n2);

        int** a21b12 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a21b12,a1+n2,a2+0,b1+0,b2+n2,n2);

        int** a22b22 = createMtrx(n2, n2);
        mtrxMulRecur(a,b,a22b22,a1+n2,a2+n2,b1+n2,b2+n2,n2);

        int** r11 = createMtrx(n2, n2);
        addMtrx(a11b11, a12b21, r11, n2);
        int** r12 = createMtrx(n2, n2);
        addMtrx(a11b12, a12b22, r12, n2);
        int** r21 = createMtrx(n2, n2);
        addMtrx(a21b11, a22b21, r21, n2);
        int** r22 = createMtrx(n2, n2);
        addMtrx(a21b12, a22b22, r22, n2);

        for(int i=0;i<n2;i++){
            for(int j=0;j<n2;j++){
                r[i][j] = r11[i][j];
            }
        }
        for(int i=0;i<n2;i++){
            for(int j=0;j<n2;j++){
                r[i][n2 + j] = r12[i][j];
            }
        }
        for(int i=0;i<n2;i++){
            for(int j=0;j<n2;j++){
                r[n2 + i][j] = r21[i][j];
            }
        }
        for(int i=0;i<n2;i++){
            for(int j=0;j<n2;j++){
                r[n2 + i][n2 + j] = r22[i][j];
            }
        }
        deleteMtrx(a11b11, n2);
        deleteMtrx(a12b21, n2);
        deleteMtrx(a11b12, n2);
        deleteMtrx(a12b22, n2);
        deleteMtrx(a21b11, n2);
        deleteMtrx(a22b21, n2);
        deleteMtrx(a21b12, n2);
        deleteMtrx(a22b22, n2);

        deleteMtrx(r11, n2);
        deleteMtrx(r12, n2);
        deleteMtrx(r21, n2);
        deleteMtrx(r22, n2);

        return;
    }
    
    // If n == 2
    r[0][0] = a[a1][a2]*b[b1][b2] + a[a1][a2 + 1]*b[b1 + 1][b2];
    r[0][1] = a[a1][a2]*b[b1][b2 + 1] + a[a1][a2 + 1]*b[b1 + 1][b2 + 1];
    r[1][0] = a[a1 + 1][a2]*b[b1][b2] + a[a1 + 1][a2 + 1]*b[b1 + 1][b2];
    r[1][1] = a[a1 + 1][a2]*b[b1][b2 + 1] + a[a1 + 1][a2 + 1]*b[b1 + 1][b2 + 1];
    //prnMtrx(r, n, n);
    return;
}

void strassens(int** a, int** b, int** result, int a1, int a2, int b1, int b2, int n){
    if(n == 2){
        int p[7];
        // p[0] = a[a1][a2] * (b[b1][b2+1] - b[b1 + 1][b2+1]);
        // p[1] = (a[a1][a2] + a[a1+0][a2+1])*b[b1+1][b2+1];
        // p[2] = (a[a1+1][a2] + a[a1 +1][a2+1])*(b[b1+0][b2+0]);
        // p[3] = (a[a1+1][a2+1])*(b[b1+1][b2+0] - b[b1+0][b2+0]); //
        // p[4] = (a[a1+0][a2+0] + a[a1+1][a2+1])*(b[b1+0][b2+0] + b[b1+1][b2+1]);//
        // p[5] = (a[a1+0][a2+1] - a[a1+1][a2+1])*(b[b1+1][b2+0] + b[b1+1][b2+1]);
        // p[6] = (a[a1+0][a2+0] - a[a1+1][a2+0])*(b[b1+0][b2+0] + b[b1+0][b2+1]);

        p[0] = (a[a1][a2] + a[a1+1][a2+1])*(b[b1][b2] + b[b1+1][b2+1]);
        p[1] = (a[a1+1][a2] + a[a1+1][a2+1])* b[b1][b2];
        p[2] = a[a1][a2] * (b[b1][b2+1] - b[b1+1][b2+1]);
        p[3] = a[a1+1][a2 + 1] * (b[b1+1][b2] - b[b1][b2]);
        p[4] = (a[a1][a2] + a[a1][a2+1])*b[b1+1][b2+1];
        p[5] = (a[a1+1][a2] - a[a1][a2])*(b[b1][b2] + b[b1][b2+1]);
        p[6] = (a[a1][a2+1] - a[a1+1][a2+1])*(b[b1+1][b1] + b[b1+1][b2+1]);


        result[0][0] = p[0] + p[3] -p[4] +p[6];//p[4] + p[3] - p[1] + p[5];
        result[0][1] = p[2] + p[4];//p[0] + p[1];
        result[1][0] = p[1] + p[3];//p[2] + p[3];
        result[1][1] = p[0] + p[2] - p[1] + p[5];//p[0] + p[4] - p[2] + p[6];
        return;
    }
    int n2 = n/2;

    int** p1 = createMatrix(n2);
    int** p2 = createMatrix(n2);
    int** p3 = createMatrix(n2);
    int** p4 = createMatrix(n2);
    int** p5 = createMatrix(n2);
    int** p6 = createMatrix(n2);
    int** p7 = createMatrix(n2);

    // P1 Calc
    int** a11a22 = createMatrix(n2);
    addMatrices(a, a, a11a22, a1+0, a2+0, a1+n2, a2+n2, n2,0);  
    int** b11b22 = createMatrix(n2);
    addMatrices(b, b, b11b22, b1+0, b2+0, b1+n2, b2+n2, n2,0);
    strassens(a11a22, b11b22, p1, 0,0, 0, 0, n2);

    //P2 Calc
    int** a21a22 = createMatrix(n2);
    addMatrices(a, a, a21a22, a1+n2, a2+0, a1+n2, a2+n2,n2, 0);
    strassens(a21a22, b, p2,0,0, b1+0,b2+0,n2);

    //P3 Calc
    int** b12b22 = createMatrix(n2);
    addMatrices(b, b, b12b22, b1+0, b2+n2, b1+n2, b2+n2,n2, 1);
    strassens(a, b12b22, p3,a1+0,a2+0,0,0,n2);

    //P4 Calc
    int** b21b11 = createMatrix(n2);
    addMatrices(b,b ,b21b11, b1+n2, b2+0, b1+0,b2+0,n2,1);
    strassens(a, b21b11, p4,a1+n2,a2+n2, 0,0, n2);

    //P5 Calc
    int** a11a12 = createMatrix(n2);
    addMatrices(a,a, a11a12,a1+0,a2+0,a1+0,a2+n2, n2, 0);
    strassens(a11a12, b, p5,0,0,b1+n2,b2+n2,n2);

    //P6 Calc
    int** a21a11 = createMatrix(n2);
    addMatrices(a,a, a21a11,a1+n2,a2+0,a1+0,a2+0,n2,1);
    int** b11b12 = createMatrix(n2);
    addMatrices(b,b, b11b12,b1+0,b2+0,b1+0,b2+n2,n2,0);
    strassens(a21a11, b11b12, p6, 0,0,0,0,n2);

    //P7 Calc
    int** a12a22 = createMatrix(n2);
    addMatrices(a,a, a12a22,a1+0, a2+n2, a1+n2,a2+n2,n2,1);
    int** b21b22 = createMatrix(n2);
    addMatrices(b,b, b21b22,b1+n2, b2+0, b1+n2,b2+n2,n2,0);
    strassens(a12a22, b21b22, p7,0,0,0,0, n2);

    int** C11 = createMatrix(n2);
    int** P1P4 = createMatrix(n2);
    int** P5P7 = createMatrix(n2);
    addMatrices(p1, p4, P1P4,0,0,0,0,n2, 0);
    addMatrices(p7, p5, P5P7,0,0,0,0,n2, 1);
    addMatrices(P1P4, P5P7, C11, 0,0,0,0, n2, 0);

    int** C12 = createMatrix(n2);
    addMatrices(p3, p5, C12,0,0,0,0, n2, 0);

    int** C21 = createMatrix(n2);
    addMatrices(p2, p4, C21,0,0,0,0, n2, 0);

    int** C22 = createMatrix(n2);
    int** P1P3 = createMatrix(n2);
    int** P2P6 = createMatrix(n2);
    addMatrices(p1, p3, P1P3,0,0,0,0,n2, 0);
    addMatrices(p6, p2, P2P6,0,0,0,0,n2, 1);
    addMatrices(P1P3, P2P6, C22, 0,0,0,0, n2, 0);

    for(int i=0;i<n2;i++){
        for(int j=0;j<n2;j++){
            result[i][j] = C11[i][j];
        }
    }

    for (int i = 0; i < n2; ++i)
    {
        for(int j=0;j < n2;j++){
            result[i][n2 + j] = C12[i][j];
        }
    }
    for(int i=0;i<n2;i++){
        for(int j=0;j<n2;j++){
            result[n2 + i][j] = C21[i][j];
        }
    }
    for(int i=0;i<n2;i++){
        for(int j=0;j<n2;j++){
            result[n2+i][n2+j] = C22[i][j];
        }
    }

    deleteMtrx(p1, n2);
    deleteMtrx(p2, n2);
    deleteMtrx(p3, n2);
    deleteMtrx(p4, n2);
    deleteMtrx(p5, n2);
    deleteMtrx(p6, n2);
    deleteMtrx(p7, n2);
    deleteMtrx(C11, n2);
    deleteMtrx(C12, n2);
    deleteMtrx(C21, n2);
    deleteMtrx(C22, n2);

    deleteMtrx(a11a22, n2);
    deleteMtrx(a21a22, n2);
    deleteMtrx(a11a12, n2);
    deleteMtrx(a21a11, n2);
    deleteMtrx(b11b22, n2);
    deleteMtrx(b12b22, n2);
    deleteMtrx(b11b12, n2);
    deleteMtrx(b21b22, n2);

    deleteMtrx(P1P4, n2);
    deleteMtrx(P5P7, n2);
    deleteMtrx(P1P3, n2);
    deleteMtrx(P2P6, n2);
    return;

}

void verify(int** rn, int** r, int n){
    int b = 0;
    int** t = createMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            t[i][j] = rn[i][j] - r[i][j];
        }
    }
    // printf("\n");
    // prnMtrx(t, n, n);
    // printf("\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(rn[i][j] != r[i][j]){
                printf("\nNot working.\n");
                return;
            }
        }
    }
    printf("Working.\n");
    deleteMtrx(t, n);
}


void test(int n){
    int** a = createMatrix(n);
    int** b = createMatrix(n);
    int** rn = createMatrix(n);
    int** rr = createMatrix(n);
    int** r = createMatrix(n);

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            a[i][j] = rand()%10;
            b[i][j] = rand()%10;
        }
    }
    printf("%d\t", n);
    clock_t begin = clock();
    // mtrxMulRecur(a,b,rn,0,0,0,0,n);
    // mtrxMul(a,b,rn,n);
    clock_t end = clock();
    double t = ((double)(end - begin))/CLOCKS_PER_SEC;
    // printf("Normal : \n");
    // prnMtrx(rn,n,n);
    // printf("Time taken = %lfs\n", t);


    // printf("\nRecursive: \n");
    begin = clock();
    mtrxMulRecur(a, b, rr, 0,0,0,0,n);
    end = clock();
    t = ((double)(end - begin))/CLOCKS_PER_SEC;
    // verify(rr, rn, n);
    // printf("Time taken = %lfs\n", t);
    printf("%lf\t", t);

    // printf("\nStrassens: \n");
    begin = clock();
    strassens(a, b, r, 0,0,0,0,n);
    end = clock();
    t = ((double)(end - begin))/CLOCKS_PER_SEC;
    // verify(rn, r, n);
    // prnMtrx(r, n,n);
    // printf("Time taken = %lfs\n", t);

    printf("%lf\n", t);
    deleteMtrx(a,n);
    deleteMtrx(b,n);
    deleteMtrx(r,n);
    deleteMtrx(rr,n);
    deleteMtrx(rn,n);
}

int main()
{
    srand(time(NULL));
    
    for(int i=1;i<=8;i++){
        test((int)pow(2,i));
    }

    return 0;
}
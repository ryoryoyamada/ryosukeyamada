#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define NUC_NUM 4
#define SIKII 5
#define RANDOMSUM 10

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

int g_hindo[NUC_NUM][MAX_SEQ_NUM];
int hindotable(int seq_num);
int hittable(int k, int gene_num, char *str, int i);
float pi[NUC_NUM][MAX_SEQ_NUM];
float si[NUC_NUM][MAX_SEQ_NUM];
float hit[MAX_SEQ_NUM][BUFSIZE];
float q[4];
float hitsum = 0.0;
float hitsumtwo = 0.0;
float num = 0.0;

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

char randomseq[RANDOMSUM][BUFSIZE];

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

int main(int argc, char* argv[]){

  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int k = hindotable(seq_num);
  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
    }
    printf("\n"); 
    for(int i = 0; i < NUC_NUM; i++){
  for(int j = 0; j < k; j++){
    printf("{%d}", g_hindo[i][j]);
  }
  printf("\n");
    }
  for(int i = 0; i < NUC_NUM; i++){
  for(int j = 0; j < k; j++){
    printf("{%f}", si[i][j]);
  }
  printf("\n");
    }
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  
  for(int i = 0; i < gene_num; i++){
  int l = hittable(k, gene_num, g_pro[i].seq, i);
  for(int j = 0; j < l; j++){
    if(hit[i][j] >= SIKII){
        printf("pro:%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
        printf("pos:%d\n", j+1);
        printf("hit(");
      for(int a = j; a < j+k; a++){
        printf("%c", g_pro[i].seq[a]);
    }
    printf(")=%.2f", hit[i][j]);
    printf("\n\n");
  }
  hit[i][j] = 0;
    }
  }
  srand((unsigned int)time(NULL));
  for(int i=0; i<10; i++){
  for(int j=0; j<400; j++){
 float x = (float)rand() / (float)RAND_MAX;
if(x <q[0]){
  randomseq[i][j] = 'A';
}
else if(x >= q[0] && x < (q[0] + q[1])){
  randomseq[i][j] = 'C';
}
else if(x >= (q[0] + q[1]) && x < (q[0] + q[1] + q[2])){
  randomseq[i][j] = 'G';
}
else if(x >= (q[0] + q[1] + q[2])){
  randomseq[i][j] = 'T';
}
}
}
printf("\n");


  for(int i = 0; i < 10; i++){
  int l = hittable(k, gene_num, randomseq[i], i);
  for(int j = 0; j < l; j++){
     hitsum += hit[i][j];
     hitsumtwo += hit[i][j] * hit[i][j];
     num++;
  }
  }
  float average = hitsum / num; 
  float std = sqrt(hitsumtwo / num - average * average);
  printf("average:%f\n", average);
  printf("std:%f\n", std);



  return 0;
}



int hindotable(int seq_num)
{
    int a = 7519429;
    int b = 4637676;
    int c = 4637676;
    int d = 7519429;
    float total = a + b + c + d;
    q[0] = a / total;
    q[1] = b / total;
    q[2] = c / total;
    q[3] = d / total;
    int k = 0;
for(int i = 0; i < seq_num; i++){
for(int j = 0; j < BUFSIZE; j++){
    if(g_motif[i][j] == 'A'){
        g_hindo[0][j]++;
    }
    else if(g_motif[i][j] == 'C'){
        g_hindo[1][j]++;
    }
    else if(g_motif[i][j] == 'G'){
        g_hindo[2][j]++;
    }
    else if(g_motif[i][j] == 'T'){
        g_hindo[3][j]++;
    }
    else{
        k = j;
        break;
    }
  } 
  }
  float tate = seq_num + NUC_NUM;
  float plus[NUC_NUM][k];
  for(int i = 0; i < NUC_NUM; i++){
  for(int j = 0; j < k; j++){
    plus[i][j] = g_hindo[i][j] + 1;
    pi[i][j] = plus[i][j] / tate;
    si[i][j] = log(pi[i][j]/q[i]);
  } 
  }
return k;
}

int hittable(int k, int gene_num, char *str, int i)
{
int l = 0;
for(l = 0; l <= strlen(str) - k; l++){
  int m = 0;
for(int j = l; j < l + k; j++){
    if(str[j] == 'A'){
        hit[i][l] += si[0][m];
    }
    else if(str[j] == 'C'){
        hit[i][l] += si[1][m];
    }
    else if(str[j] == 'G'){
        hit[i][l] += si[2][m];
    }
    else if(str[j] == 'T'){
        hit[i][l] += si[3][m];
    }
    m++;
  } 
}
  return l;
}


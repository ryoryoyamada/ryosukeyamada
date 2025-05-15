
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

int g_hindo[4][20];
int hindotable(int seq_num);
float pi[4][20];
float si[4][20];
float hit[BUFSIZE];

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

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
  int l = hindotable(seq_num);
  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
    }
    printf("\n"); 
    for(int j = 0; j < 4; j++){
  for(int k = 0; k < l; k++){
    printf("{%d}", g_hindo[j][k]);
  }
  printf("\n");
    }
    for(int j = 0; j < 4; j++){
  for(int k = 0; k < l; k++){
    printf("{%f}", si[j][k]);
  }
  printf("\n");
    }
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  int m = hittable(l, gene_num);
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
    printf("%f\n", hit[i]);
  }
  
  return 0;
}

int hindotable(int seq_num)
{
    int a = 7519429;
    int b = 4637676;
    int c = 4637676;
    int d = 7519429;
    float total = a + b + c + d;
    float q[4];
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
  float tate = seq_num + 4;
  float plus[4][k];
  for(int i = 0; i < 4; i++){
  for(int j = 0; j < k; j++){
    plus[i][j] = g_hindo[i][j] + 1;
    pi[i][j] = plus[i][j] / tate;
    si[i][j] = log(pi[i][j]/q[i]);
  } 
  }
return k;
}

int hittable(int k, int gene_num)
{
int l = 0;
for(int i = 0; i < gene_num; i++){
for(int j = l; j < l + k; j++){
    if(g_pro[i].seq[j] == 'A'){
        hit[l] = + si[0][j];
    }
    else if(g_pro[i].seq[j] == 'C'){
        hit[l] = + si[1][j];
    }
    else if(g_pro[i].seq[j] == 'G'){
        hit[l] = + si[2][j];
    }
    else if(g_pro[i].seq[j] == 'T'){
        hit[l] = + si[3][j];
    }
    if(l >= gene_num - k){
        break;
    }
    l++;
  } 
}
  return l;
}
/*
 * inputRef.cpp
 *
 *  Created on: 2021年5月18日
 *      Author: lin
 */

#include "inputRef.h"

void ReadSeq(char **seq1, uint64_t *seq_length, char *p_ref) {
	uint32_t buffer_size = 256;
	char buffer_line[256];
	//初始化
	memset(buffer_line, 0, buffer_size);

	FILE *fp;
	fp = fopen(p_ref, "r+");
	if (fp == NULL) {
		cout << "file can not be open!" << endl;
		return;
	}

	uint64_t total_size = 0;
	//指针定位到文件尾
	fseek(fp, 0, 2);
	//计算序列总长度
	total_size = ftell(fp);

	char *seq;
	seq = (char*) malloc(sizeof(char) * total_size);

	//指针定位到文件头
	fseek(fp, 0, 0);

	uint64_t len = 0;
	while (fgets(buffer_line, buffer_size - 1, fp) != NULL) {
		if (buffer_line[0] == '>')
			continue;
		else {
			for (uint32_t i = 0; i < buffer_size; i++) {
				if (buffer_line[i] == '\n' || buffer_line[i] == '\0') {
					break;
				}
				if (buffer_line[i] >= 'a') {
					buffer_line[i] -= 32;
				}
				seq[len] = buffer_line[i];
				len++;
			}
		}
		memset(buffer_line, 0, buffer_size);
	}
	*seq_length = len;
	*seq1 = seq;
	cout << "the length of seq is: " << len << endl;

}

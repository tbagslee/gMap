/*
 * Main.cpp
 *
 *  Created on: 2021年5月17日
 *      Author: lin
 */
#include "basic.h"
#include "inputReads.h"
#include "inputRef.h"
#include "seeding.h"
#include "mappers.h"
#include"outputSam.h"
uint32_t num_read;

int main(int argc, char **argv) {

	//当前可执行文件下文件名read
	char *p_read_path;

	//ref参考基因
	char *p_ref_path;

//	索引文件
	char *index_file_path;

	//输出sam
	char *output_file_path;

	//read的长度
	uint32_t len_read;

	//任务量（读取几个read）
	uint32_t task_size;

	//阈值,默认为1
	uint32_t e = 1;

	//线程
	uint32_t thread_num;

	//读取命令行参数
	for (int32_t i = 1; i < argc; i = i + 2) {
		if (argv[i][0] == '-' && argv[i][1] == 'r') {
			p_read_path = argv[i + 1];
		}
		if (argv[i][0] == '-' && argv[i][1] == 'R') {
			p_ref_path = argv[i + 1];
		}
		if (argv[i][0] == '-' && argv[i][1] == 'L') {
			len_read = atoi(argv[i + 1]);
		}
		if (argv[i][0] == '-' && argv[i][1] == 'T') {
			task_size = atoi(argv[i + 1]);
		}
		if (argv[i][0] == '-' && argv[i][1] == 'e') {
			e = atoi(argv[i + 1]);
		}
		if (argv[i][0] == '-' && argv[i][1] == 't') //
				{
			thread_num = atoi(argv[i + 1]);
		}
		if (argv[i][0] == '-' && argv[i][1] == 'o') {
			output_file_path = argv[i + 1];
		}
//		if(argv[i][0] == '-' && argv[i][1] == 'i') {
//			index_file_path = argv[i + 1];
//		}
	}

	//开始时间 结束时间
	struct timeval tvs, tve;
	gettimeofday(&tvs, NULL);
	cout << "start..." << endl;

	//调用input_method方法的例子
	char *p_reads;
	uint32_t start_line = 0;
	input_reads(&p_reads, &num_read, p_read_path, start_line, task_size,
			len_read);
//	cout << p_reads <<endl;
//	cout << num_read <<endl;

//調用ref的例子
	char *p_ref;
	uint64_t ref_length;
	ReadSeq(&p_ref, &ref_length, p_ref_path);
//	cout << p_ref <<endl;
//	cout << ref_length <<endl;

//seed_indexing调用
	struct NodeBit **p_BplusTrees;
	vector<struct seed_message> VS_message;
	uint32_t *seed_find_read;
	p_BplusTrees = seed_indexing(p_reads, num_read, len_read, thread_num, e,
			VS_message, &seed_find_read);
	cout << VS_message.size() << endl;
//	for (uint32_t i = 0; i < VS_message.size(); i++) {
//		//输出种子类型的个数p_read_length
//		cout << VS_message[i].seed_type_count << endl;
//	}
//	cout << "test_saInterval" << endl;
//	for (uint32_t i = 0; i < VS_message.size(); i++) {
//		//输出SA
//		cout << VS_message[i].saInterval[0] << endl;
//		cout << VS_message[i].saInterval[1] << endl;
//	}
//
//	cout << "seed_sa" << endl;
//	for (uint32_t i = 0; i < VS_message.size(); i++) {
//		for (uint32_t k = 0;k < VS_message[i].saInterval[1] - VS_message[i].saInterval[0] + 1; k++) {
//			cout << VS_message[i].seed_sa[k] << endl;
//		}
//	}
//
//	cout << "test_seed_find_read" <<endl;
//	for(uint32_t i = 0;i < num_read;i++){
//		for(uint32_t j = 0;j < e+1;j++){
//			cout << *(seed_find_read + i * (e + 1) + j) << endl;
//		}
//	}

	//检测message中候选区间size
	cout << "size:" << VS_message.size() <<endl;
	uint32_t num = 0;
	for(uint32_t i = 0;i<VS_message.size();i++){
		if(VS_message[i].saInterval[1]>=VS_message[i].saInterval[0]){
			num+=VS_message[i].saInterval[1]-VS_message[i].saInterval[0]+1;
		}
	}
	cout << "num:" << num <<endl;

	struct naiveMapper_para p;
	p.VS_message=VS_message;
	p.e=e;
	p.len_read=len_read;
	p.p_BplusTrees=p_BplusTrees;
	p.p_reads=p_reads;
	p.p_ref=p_ref;
	p.ref_length=ref_length;
	p.task_size=task_size;
	p.thread_num=thread_num;
	p.seed_find_read=seed_find_read;
	nativeMapper(p);

	//    // 输出 Sam

	OutSam *output_sam;
	char* read_qual = NULL;
	outSam(output_sam,output_file_path, p_ref_path,e, \
							0, 0, p_read_path, strlen(p_read_path),\
							p_reads,read_qual,len_read);

//	struct OutputQueue *output_queue;
//	  // Load reference
//	SequenceBatch reference_sequence_batch;
//    initialize_sequence_batch(&reference_sequence_batch);
//    initialize_sequence_batch_loading(p_ref_path, &reference_sequence_batch);
//    load_all_sequences_into_sequence_batch(&reference_sequence_batch);
//   // initialize output_queue
//    uint32_t output_queue_max_size = 100000;
//    initialize_output_queue(output_file_path, &reference_sequence_batch, thread_num, output_queue_max_size, output_queue);
////	void *p;  //指向指针函数
//	void *output_queue_v;
//	output_queue_v = (void*)malloc(sizeof(struct OutputQueue));
//	output_queue_v = &output_queue;
//	fprintf(stderr, "output: %s\n", output_file_path);
//	output_sam(output_queue);

	cout << "end..." << endl;
	gettimeofday(&tve, NULL);
	double span = tve.tv_sec - tvs.tv_sec
			+ (tve.tv_usec - tvs.tv_usec) / 1000000.0;
	cout << "seeding time is: " << span << endl;
}

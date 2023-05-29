/*
 * seeding.cpp
 *
 *  Created on: 2021年5月18日
 *      Author: lin
 */

#include "seeding.h"

void* Parallel_cal_seeds(void *p) {
//	cout << "test-cal" << endl;
	struct bit256KmerPara para_tmp;

	struct para_multi_seed *tmp = (struct para_multi_seed*) p;
	uint64_t *seed_code;
	struct seed_message *seed_message;
	uint8_t *seed_thread_id;
	uint32_t seed_num = 1; // the seed_num value is 1 bigger than the actual size of the seeds
	uint32_t seed_64num = tmp->len_read / (tmp->e + 1) / 32 + 1 + 1;
	seed_code = (uint64_t*) malloc(sizeof(uint64_t) * seed_num * seed_64num);
	seed_message = (struct seed_message*) malloc(
			sizeof(struct seed_message) * seed_num);
	seed_thread_id = (uint8_t*) malloc(sizeof(uint8_t) * seed_num);

	//从哪个read开始和结束
	for (uint32_t i = tmp->start; i <= tmp->end; i++) {
		struct seed_message message_tmp;
		uint32_t *p_SAI;
		char seed_tmp[64];
		//i*len_read表示当前read的位置
		for (uint32_t j = 0; j <= tmp->e; j++) {
			//1)划分
			uint32_t start_pos;
			uint32_t seed_len;
			generate_seed(&start_pos, &seed_len, tmp->len_read, tmp->e, j);
			get_para(&para_tmp, seed_len);

			uint32_t start_tmp;
			start_tmp = (seed_num - 1) * seed_64num;
			seed_code[start_tmp] = seed_len;
			cal_hash_value_directly_256bit(
					tmp->p_reads + i * tmp->len_read + start_pos,
					seed_code + start_tmp + 1, para_tmp);
			seed_code = (uint64_t*) realloc(seed_code,
					sizeof(uint64_t) * (seed_num + 1) * seed_64num);

			message_tmp.seed_len = seed_len;
			message_tmp.e = tmp->e;
			//seed_type_count暂时存read的id，后面sav方法可以给seed_message的read_id
			message_tmp.seed_type_count = i;
			//message_tmp的p_read_id和seed_type_count还未算出来

			strncpy(seed_tmp, tmp->p_reads + i * tmp->len_read + start_pos,
					seed_len);
			p_SAI = calc_SArangeSeq(tmp->FMidx, seed_tmp);
			message_tmp.saInterval[0] = p_SAI[0];
			message_tmp.saInterval[1] = p_SAI[1];

			seed_message[seed_num - 1] = message_tmp;
			seed_message = (struct seed_message*) realloc(seed_message,
					sizeof(struct seed_message) * (seed_num + 1));

			//后面sav方法做判断用，对线程id进行加密？
			seed_thread_id[seed_num - 1] = (bit256hashFFunction(
					seed_code + start_tmp + 1, para_tmp)) % tmp->thread_num;
			seed_thread_id = (uint8_t*) realloc(seed_thread_id,
					sizeof(uint8_t) * (seed_num + 1));
//			cout << seed_thread_id[seed_num - 1] << endl;

			seed_num = seed_num + 1;
		}
	}
	*(tmp->p_seed) = seed_code;
	*(tmp->p_seed_message) = seed_message;
	*(tmp->p_seed_thread_id) = seed_thread_id;
	*(tmp->p_seed_num) = seed_num - 1;
}

void* Parallel_sav_seeds(void *p) {
//	cout << "test-sav" << endl;
	struct bit256KmerPara para_tmp;

	struct para_multi_seed_final *tmp = (struct para_multi_seed_final*) p;
	uint32_t num_id = 0;
	for (uint32_t i = 0; i < tmp->thread_num; i++) {
		for (uint32_t j = 0; j < tmp->p_seed_num[i]; j++) {
			if (tmp->seed_thread_id[i][j] == tmp->thread_id) {
				//计算每个种子的哈希
				get_para(&para_tmp, tmp->seed_message[i][j].seed_len);
				para_tmp.kmer64Len += 1;
				struct nodeBit c_tmp_hashtable;
				c_tmp_hashtable.hashValue = tmp->p_seed[i]
						+ j * para_tmp.kmer64Len;
				int64_t x;
				x = getHashFTableValue(tmp->p_tmp, c_tmp_hashtable.hashValue,
						para_tmp);
				if (x != -1) {
					(*(tmp->VS_message))[x].p_read_id = (uint32_t*) realloc(
							(*(tmp->VS_message))[x].p_read_id,
							sizeof(uint32_t)
									* ((*(tmp->VS_message))[x].seed_type_count
											+ 1));
					(*(tmp->VS_message))[x].p_read_id[(*(tmp->VS_message))[x].seed_type_count] =
							tmp->seed_message[i][j].seed_type_count;
					(*(tmp->VS_message))[x].seed_type_count++;
					*(tmp->seed_find_read + num_id + j) = x;
				} else {
					struct seed_message tmp_message;
					tmp_message.e = tmp->seed_message[i][j].e;
					tmp_message.seed_len = tmp->seed_message[i][j].seed_len;
					tmp_message.p_read_id = (uint32_t*) malloc(
							sizeof(uint32_t));
					tmp_message.seed_type_count = 0;
					//上面的cal方法暂时将read的id赋给了seed_type_count
					tmp_message.p_read_id[tmp_message.seed_type_count] =
							tmp->seed_message[i][j].seed_type_count;
					tmp_message.seed_type_count++;

					tmp_message.saInterval[0] =
							tmp->seed_message[i][j].saInterval[0];
					tmp_message.saInterval[1] =
							tmp->seed_message[i][j].saInterval[1];

					//calc_SA
					if (tmp_message.saInterval[0] <= tmp_message.saInterval[1]) {
						tmp_message.seed_sa =(uint32_t*) malloc(sizeof(uint32_t) * (tmp_message.saInterval[1] - tmp_message.saInterval[0]
														+ 1));
						for (uint32_t m = 0;m < tmp_message.saInterval[1] - tmp_message.saInterval[0] + 1; m++) {
							tmp_message.seed_sa[m] = calc_SA(tmp->FMidx,tmp_message.saInterval[0] + m);
						}
					}

					pthread_mutex_lock(tmp->m);
					//与单线程相比，种子放入B+树的时候生成的id号可能不同，但对应位置的种子是完全一样的
					*(tmp->seed_find_read + num_id + j) = *(tmp->array_id);
					c_tmp_hashtable.arrayID = *(tmp->array_id);
					*(tmp->array_id) = *(tmp->array_id) + 1;
					(*(tmp->VS_message)).push_back(tmp_message);
					pthread_mutex_unlock(tmp->m);
					bit256insertHashFTable(tmp->p_tmp, c_tmp_hashtable,
							para_tmp);
				}
			}
		}
		num_id += tmp->p_seed_num[i];
	}
}

void generate_seed(uint32_t *start_pos, uint32_t *seed_len, uint32_t len_read,
		uint32_t e, uint32_t seed_id) {
	//read长度100,划分种子为34,33,33
	uint32_t p = len_read / (e + 1);
	uint32_t q = len_read % (e + 1);
	if (seed_id < q) {
		*seed_len = p + 1;
		*start_pos = (p + 1) * seed_id;
	} else {
		*seed_len = p;
		*start_pos = (p + 1) * q + p * (seed_id - q);
	}

}

struct NodeBit** seed_indexing(char *p_reads, uint32_t num_read,
		uint32_t len_read, uint32_t thread_num, uint32_t e,
		vector<struct seed_message> &VS_message, uint32_t **seed_find_read) {

	//定义seed找到read的数组大小
	uint32_t *s_f_r;
	s_f_r = (uint32_t*) malloc(sizeof(uint32_t) * num_read * (e + 1));

	//FM-index计算
	sFMindex FMidx;
	char *path = (char*) malloc(sizeof(char) * 20);
	path = ".";
	read_bfile2index(path, FMidx, 0);

	struct bit256KmerPara para_tmp;

	//(1) 初始化B+树
	struct NodeBit **p_tmp;
	p_tmp = bit256initialHashFTable();

	//(2) 产生种子并将它们放在B+树
	//单线程
	if (thread_num == 1) {

		//B+树上的ID，不重复
		uint64_t arrayId = 0;

		//每次处理一个read
		for (uint32_t i = 0; i < num_read; i++) {

			//放入每個seed的长度和哈希
			uint64_t *seed_code;
			seed_code = (uint64_t*) malloc(
					sizeof(uint64_t) * (len_read / (e + 1) / 32 + 1 + 1));

			//存储位置
			char seed_tmp[64];
			uint32_t *p_SAI;
			//每次处理read上的一个seed
			for (uint32_t j = 0; j <= e; j++) {
				//划分
				uint32_t start_pos;
				uint32_t seed_len;
				generate_seed(&start_pos, &seed_len, len_read, e, j);

				//对seedcode进行赋值
				get_para(&para_tmp, seed_len);
				seed_code[0] = seed_len;
				cal_hash_value_directly_256bit(
						p_reads + i * len_read + start_pos, seed_code + 1,
						para_tmp);
				para_tmp.kmer64Len += 1;

				//输出每个种子和每个种子的长度和哈希-测试用
//				for (uint32_t i_test = i * len_read + start_pos;
//						i_test < i * len_read + start_pos + seed_len;
//						i_test++) {
//					cout << p_reads[i_test];
//				}
//				cout << endl;
//				cout << seed_code[0] << endl;
//				cout << seed_code[1] << endl;

				//插入
				struct nodeBit c_tmp_hashtable;
				c_tmp_hashtable.hashValue = seed_code;
				int64_t x;
				x = getHashFTableValue(p_tmp, c_tmp_hashtable.hashValue,
						para_tmp);
				if (x != -1) {
					//如果当前种子存在于B+树上
					VS_message[x].p_read_id = (uint32_t*) realloc(
							VS_message[x].p_read_id,
							sizeof(uint32_t)
									* (VS_message[x].seed_type_count + 1));
					VS_message[x].p_read_id[VS_message[x].seed_type_count] = i;
					VS_message[x].seed_type_count++;
					*(s_f_r + i * (e + 1) + j) = x;
				} else {

					//如果當前种子不在B+树上
					c_tmp_hashtable.arrayID = arrayId;
					//保存当前种子的信息
					strncpy(seed_tmp, p_reads + i * len_read + start_pos,
							seed_len);
					p_SAI = calc_SArangeSeq(FMidx, seed_tmp);

					struct seed_message tmp;
					tmp.e = 0;
					tmp.seed_len = seed_len;
					tmp.p_read_id = (uint32_t*) malloc(sizeof(uint32_t));
					tmp.seed_type_count = 0;
					tmp.p_read_id[tmp.seed_type_count] = i;
					tmp.seed_type_count++;
					tmp.saInterval[0] = p_SAI[0];
					tmp.saInterval[1] = p_SAI[1];
					if (tmp.saInterval[1] >= tmp.saInterval[0]) {
						tmp.seed_sa = (uint32_t*) malloc(
								sizeof(uint32_t) * (tmp.saInterval[1] - tmp.saInterval[0] + 1));
						for (uint32_t k = 0; k < tmp.saInterval[1] - tmp.saInterval[0] + 1; k++) {
							tmp.seed_sa[k] = calc_SA(FMidx, tmp.saInterval[0] + k);
						}
					}

					VS_message.push_back(tmp);
					bit256insertHashFTable(p_tmp, c_tmp_hashtable, para_tmp);
					*(s_f_r + i * (e + 1) + j) = arrayId;
					arrayId++;
				}
			}
		}
	} else {
		//多线程
		pthread_t *t;
		t = (pthread_t*) malloc(sizeof(pthread_t) * thread_num);
		pthread_mutex_t m;
		pthread_mutex_init(&m, NULL);

		//1 计算每个线程处理多少个read
		uint32_t p = num_read / thread_num;
		uint32_t q = num_read % thread_num;

		//多线程cal时带的参数
		struct para_multi_seed *para;
		para = (struct para_multi_seed*) malloc(
				sizeof(struct para_multi_seed) * thread_num);

		uint32_t cur = 0;
		uint64_t *p_seed[thread_num];
		struct seed_message *seed_message[thread_num];
		uint8_t *seed_thread_id[thread_num];
		uint32_t p_seed_num[thread_num];
		uint32_t array_id = 0;

		for (uint32_t i = 0; i < thread_num; i++) {
			para[i].FMidx = FMidx;
			para[i].start = cur;
			if (i < q) {
				cur = cur + p + 1;
			} else {
				cur = cur + p;
			}
			para[i].end = cur - 1;
			para[i].p_reads = p_reads;
			para[i].len_read = len_read;
			para[i].p_seed = &(p_seed[i]);
			para[i].p_seed_message = &(seed_message[i]);
			para[i].p_seed_thread_id = &(seed_thread_id[i]);
			para[i].p_seed_num = p_seed_num + i;
			para[i].thread_id = i;
			para[i].thread_num = thread_num;
			para[i].e = e;

			if (pthread_create(t + i, NULL, Parallel_cal_seeds,
					(void*) (para + i)) != 0) {
				cout << "error!" << endl;
			}
		}
		for (uint32_t i = 0; i < thread_num; i++) {
			pthread_join(t[i], NULL);
		}
		free(para);

		//2 插入
		struct para_multi_seed_final *para_final;
		para_final = (struct para_multi_seed_final*) malloc(
				sizeof(struct para_multi_seed_final) * thread_num);
		for (uint32_t i = 0; i < thread_num; i++) {

			para_final[i].FMidx = FMidx;
			para_final[i].VS_message = &VS_message;
			para_final[i].array_id = &array_id;
			para_final[i].m = &m;
			para_final[i].p_seed = p_seed;
			para_final[i].p_seed_num = p_seed_num;
			para_final[i].p_tmp = p_tmp;
			para_final[i].seed_message = seed_message;
			para_final[i].seed_thread_id = seed_thread_id;
			para_final[i].thread_id = i;
			para_final[i].thread_num = thread_num;
			para_final[i].seed_find_read = s_f_r;

			if (pthread_create(t + i, NULL, Parallel_sav_seeds,
					(void*) (para_final + i)) != 0) {
				cout << "error!" << endl;
			}
		}
		for (uint32_t i = 0; i < thread_num; i++) {
			pthread_join(t[i], NULL);
		}

		//释放资源
		free(para_final);
		free(t);
		for (uint32_t i = 0; i < thread_num; i++) {
			if (p_seed[i] != NULL) {
				free(p_seed[i]);
			}
			if (seed_message[i] != NULL) {
				free(seed_message[i]);
			}
			if (seed_thread_id[i] != NULL) {
				free(seed_thread_id[i]);
			}
		}
	}
	*seed_find_read = s_f_r;
	return p_tmp;
}

uint32_t get_seed_pos(uint32_t len_read, uint32_t e, uint32_t seed_id) {
	uint32_t p = len_read / (e + 1);
	uint32_t q = len_read % (e + 1);
	if (seed_id < q) {
		return (p + 1) * seed_id;
	} else {
		return (p + 1) * q + p * (seed_id - q);
	}
}

//int outputSam_main(int argc, char *argv[])
//{
//	destroy_output_queue(&output_queue);
//	destroy_input_queue(&input_queue);
//	destroy_index(&index);
//	finalize_sequence_batch_loading(&reference_sequence_batch);
//	destory_sequence_batch(&reference_sequence_batch);
//	return 0;
//}

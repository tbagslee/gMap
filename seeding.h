/*
 * seeding.h
 *
 *  Created on: 2021年5月19日
 *      Author: lin
 */

#ifndef SEEDING_H_
#define SEEDING_H_
#include"basic.h"
#include"BplusTreeBit.h"
#include"Hash.h"
#include"exactMatchFMindex.h"
#include"outputSam.h"

struct seed_message
{
	uint8_t e;
	uint8_t seed_len;
	uint32_t *p_read_id;
	uint32_t seed_type_count;
	uint32_t saInterval[2];
	uint32_t *seed_sa;
};

struct para_multi_seed
{
	sFMindex FMidx;
	uint32_t start;
	uint32_t end;
	char* p_reads;
	uint32_t len_read;
	uint64_t ** p_seed;
	uint32_t * p_seed_num;
	uint8_t ** p_seed_thread_id;
	struct seed_message **p_seed_message;
	uint32_t thread_id;
	uint32_t thread_num;
	uint32_t e;
};

struct para_multi_seed_final
{
	sFMindex FMidx;
	uint64_t ** p_seed;
	struct seed_message **seed_message;
	uint8_t **seed_thread_id;
	uint32_t * p_seed_num;
	uint32_t thread_num;
	uint32_t thread_id;
	uint32_t e;
	uint32_t *seed_find_read;
	struct NodeBit** p_tmp;
	vector <struct seed_message> *VS_message;
	pthread_mutex_t * m;
	uint32_t* array_id;
};

void generate_seed(uint32_t *start_pos, uint32_t *seed_len, uint32_t len_read,
		uint32_t e, uint32_t seed_id);

struct NodeBit** seed_indexing(char *p_reads, uint32_t num_read,
		uint32_t len_read, uint32_t thread_num, uint32_t e,
		vector<struct seed_message> &VS_message, uint32_t **seed_find_read);

uint32_t get_seed_pos(uint32_t len_read, uint32_t e, uint32_t seed_id);


#endif /* SEEDING_H_ */

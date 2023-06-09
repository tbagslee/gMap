/*
 * mapping.h
 *
 *  Created on: 2021年5月21日
 *      Author: lin
 */

#ifndef MAPPERS_H_
#define MAPPERS_H_

#include "basic.h"
#include "seeding.h"
#include "prealignment.h"

struct naiveMapper_para
{
	char *p_reads;
	uint32_t len_read;
	uint32_t task_size;
	uint32_t *seed_find_read;
	vector <struct seed_message> VS_message;

	char* p_ref;
	uint64_t ref_length;

	struct NodeBit** p_BplusTrees;

	uint32_t e;
	uint32_t thread_num;
};

struct para_multi_p{
	uint32_t start;
	uint32_t end;
	struct naiveMapper_para *p;
	pthread_mutex_t *m;
};

void* Parall_nativeMapper(void *p);

void nativeMapper(struct naiveMapper_para p);



#endif /* MAPPERS_H_ */

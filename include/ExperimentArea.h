#ifndef EXPERIMENT_AREA_H
#define EXPERIMENT_AREA_H

#include <vector>
#include <deque>
#include <string.h>
#include "GlobalDefinitions.h"

namespace ParameterPile {

	class area_vector
	{
	protected:
		bool _is_valid;
		STD_CONT<int> _vec;
		int _last_returned_index;
		STD_CONT<int>::iterator _last_returned_index_left;
	public:
		area_vector(void);
		//area_vector(STD_CONT<int> channels); //TODO:
		int get_order_index_by_index(int ind) const;
		int get_index_by_order_index(int ind) const;
		//TODO: rework in thread-save manner. (via size, get[] or altogether as iterator)
		//TODO: actually mb all settings objects should be copied for each thread?
		int get_next_index(void); //for running through all indices
		int get_next_index(int after);
		void push_pair(int left, int right);
		void push_back(int val);
		void push (int left, int right);
		void push (int val);
		bool contains(int index) const;
		bool empty(void) const;
		bool isValid(void) const;
		int &back(void);
		int &front(void);
		STD_CONT<area_vector> split_area(int N) const;
		area_vector intersect(area_vector with);
		void reset(); //clears _last_returned_index etc.
		void erase(); //clears vector
		//void refine (void); //e.g. [2,3][3,4] to [2,4] |OR| [4,5] [1,7] to [1,7]
		std::size_t real_size(void) const;
		//bool operator == (const area_vector& with) const;
	};

	class experiment_area //done //TODO - make analysis via this class. //->NextFile?
	{
	public:
		enum Type { Area, Point };
	protected:
		Type _type;
	public:
		experiment_area(Type type = Type::Area);
		experiment_area to_point(void);

		STD_CONT<std::string> experiments;
		area_vector runs; //contains pairs [from, to]
		area_vector channels; //contains pairs [from, to]
		area_vector sub_runs; //contains pairs [from, to]

		bool isValid(void);
		bool contains(ParameterPile::experiment_area what);//new draw_required
		std::size_t real_size(void);
	};
};

#endif

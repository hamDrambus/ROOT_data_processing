#ifndef EXPERIMENT_AREA_H
#define EXPERIMENT_AREA_H

#include <vector>
#include <deque>
#include <string.h>
#include "GlobalDefinitions.h"

namespace ParameterPile {

	//TODO: rework - create something like channel/analysis manifest:
	//experiment-> folder
	//          -> accepted_events.txt (with adjusted trigger)
	//          -> output_folder
	//          -> {ch1, ch2, ..., ch32}
	//             |- find_peaks
	//             |- baseline restoration times
	//             |- double_I times
	//             |- invert
	//             |- device: {"PMT", "MPPC", "GEM"}
	//             |- threshold, use_average?, calculate integral?, calculate double integral?, their parameters as well
	//             ...and so on. This way each channel may be fine-tuned. (May be filled using default for each experiment,
	//             and there is a default for all experements as well). Setting parameters should be done in separate files (190404.cpp, 180705.cpp)
	//             so that analysis can be reproduced easily just by changing which file is loaded (compiled) similar to Post_processor.
	//             The variant with using XML will be too cumbersome
	//In addition, run and subrun numbers should also be saved with peaks, because now it is just silently assumed, that events are in succession in both
	//Data_processing and Post_procesing, but this won't be the case if only events from accepted_events.txt are processed again in Data_processing.

	template <class T>
	class indexed_info { //called channel_info in Post_processing
	protected:
		std::deque<std::pair<int, T> > _data;
	public:
		indexed_info()
		{}
		~indexed_info()
		{}
		template <class U>
		bool isSameIndices(const indexed_info<U>& b) const
		{
			if (_data.size() != b._data.size())
				return false;
			else {
				for (std::size_t ind = 0, ind_end_ = _data.size(); ind < ind_end_; ++ind)
					if (_data[ind].first != b._data[ind].first)
						return false;
			}
			return true;
		}
		bool isSameIndices(const std::deque<int>& indices) const
		{
			if (_data.size() != indices.size())
				return false;
			else {
				for (std::size_t ind = 0, ind_end_ = _data.size(); ind < ind_end_; ++ind)
					if (_data[ind].first != indices[ind])
						return false;
			}
			return true;
		}
		void push(const int& index, const T& data) //preserves channel sorting
		{
			std::size_t sz = _data.size();
			if (0 == sz) {
				_data.push_back(std::pair<int, T>(index, data));
				return;
			}
			if (index < _data.front().first) {
				_data.insert(_data.begin(), std::pair<int, T>(index, data));
				return;
			}
			if (index > _data.back().first) {
				_data.push_back(std::pair<int, T>(index, data));
				return;
			}
			std::pair<std::size_t, std::size_t> inds = get_bounds(index);
			if (inds.first == inds.second) //do not insert points with equal channel, replace only
				_data[inds.first].second = data;
			else
				_data.insert(_data.begin() + inds.second, std::pair<int, T>(index, data));
		}
		void push_back(const int& index, const T& data) //ignores channel sorting
		{
			_data.push_back(std::pair<int, T>(index, data));
		}
		T* info(const int& index)
		{
			std::size_t ch_ind = get_index(index);
			if (std::numeric_limits<std::size_t>::max() == ch_ind) {
				return NULL;
			}
			return &(_data[ch_ind].second);
		}
		const T* info(const int& index) const
		{
			std::size_t ch_ind = get_index(index);
			if (std::numeric_limits<std::size_t>::max() == ch_ind) {
				return NULL;
			}
			return &(_data[ch_ind].second);
		}

		T& operator [] (const std::size_t& ch_ind)
		{
			std::size_t sz = _data.size();
			if (!(ch_ind < sz))
				throw std::out_of_range("indexed_info<T>::operator[] index is out of range");
			return _data[ch_ind].second;
		}
		void erase(const int& index)
		{
			std::size_t ch_ind = get_index(index);
			if (std::numeric_limits<std::size_t>::max() == ch_ind) {
				std::cout << "indexed_info<T>::erase: Warning! Does not contain entry " << index << std::endl;
				return;
			}
			erase_at(ch_ind);
		}
		void erase_at(const std::size_t& ch_ind)
		{
			std::size_t sz = _data.size();
			if (!(ch_ind < sz))
				throw std::out_of_range("indexed_info<T>::erase_at index is out of range");
			_data.erase(_data.begin() + ch_ind);
		}

		void clear(void)
		{
			std::deque<std::pair<int, T> >().swap(_data);
		}

		std::size_t size(void) const
		{
			return _data.size();
		}

		bool empty(void) const
		{
			return _data.empty();
		}

		bool contains(const int& index) const
		{
			if (std::numeric_limits<std::size_t>::max() != get_index(index))
				return true;
			return false;
		}

		std::size_t get_index(const int& index) const
		{
			std::size_t sz = _data.size();
			if (0 == sz) {
				return std::numeric_limits<std::size_t>::max();
			}
			if (index < _data.front().first) {
				return std::numeric_limits<std::size_t>::max();
			}
			if (index > _data.back().first) {
				return std::numeric_limits<std::size_t>::max();
			}
			std::pair<std::size_t, std::size_t> inds = get_bounds(index);
			if (inds.first == inds.second)
				return inds.first;
			else
				return std::numeric_limits<std::size_t>::max();
		}
		int index(const std::size_t& ch_ind) const
		{
			std::size_t sz = _data.size();
			if (!(ch_ind < sz))
				throw std::out_of_range("indexed_info<T>::channel index is out of range");
			return _data[ch_ind].first;
		}

		std::pair<std::size_t, std::size_t> get_bounds(const int& index) const
		{
			std::pair<std::size_t, std::size_t> out(std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max());
			std::size_t sz = _data.size();
			if (0 == sz)
				return out;
			if (index <= _data.front().first) {
				out = std::pair<std::size_t, std::size_t>(0, 0);
				return out;
			}
			if (index >= _data.back().first) {
				out = std::pair<std::size_t, std::size_t>(sz - 1, sz - 1);
				return out;
			}
			//find first ch which is not less that channel. That is index bounding X_point: chs[first] <= channel < chs[first + 1]
			//See std::lower_bound and std::upper_bound:
			//std::lower_bound(xys.begin(), xys.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b)->bool{
			//	return a.first<b.first;
			//});
			std::size_t count = sz;
			std::size_t first = 0;
			while (count > 0) {
				std::size_t step = count / 2;
				std::size_t ind = first + step;
				if (!(index < _data[ind].first)) {
					first = ++ind;
					count -= step + 1;
				} else
					count = step;
			}
			//first is such, that channel>=chs[first-1] and channel<chs[first]
			//first always != 0 here
			--first;
			if (index == _data[first].first) {
				out = std::pair<std::size_t, std::size_t>(first, first);
				return out;
			}
			out = std::pair<std::size_t, std::size_t>(first, first + 1);
			return out;
		}
		void sort(void)
		{
			std::sort(_data.begin(), _data.end(), [](const std::pair<int, T>& a, const std::pair<int, T>& b)->bool {
				return a.first < b.first;
			});
		}
		bool is_sorted(void) const
		{
			for (std::size_t ind = 0, ind_end_ = _data.size(); (ind + 1) != ind_end_; ++ind)
				if (_data[ind].first >= _data[ind + 1].first)
					return false;
			return true;
		}
	};

	template <class T>
	class accepted_events : public indexed_info<indexed_info<T> > { //stores run->subrun->adjusted trigger value
	public:
		accepted_events(void)
		{}
		~accepted_events()
		{}
		void push(const int& run, const int& subrun, const T& trigger_offset)
		{
			indexed_info<double> single_element;
			single_element.push(subrun, trigger_offset);
			std::size_t sz = this->_data.size();
			if (0 == sz) {
				this->_data.push_back(std::pair<int, indexed_info<double> >(run, single_element));
				return;
			}
			if (run < this->_data.front().first) {
				this->_data.insert(this->_data.begin(), std::pair<int, indexed_info<double> >(run, single_element));
				return;
			}
			if (run > this->_data.back().first) {
				this->_data.push_back(std::pair<int, indexed_info<double> >(run, single_element));
				return;
			}
			std::pair<std::size_t, std::size_t> inds = this->get_bounds(run);
			if (inds.first != inds.second)
				this->_data.insert(this->_data.begin() + inds.second, std::pair<int, indexed_info<double> >(run, single_element));
			else
				this->_data[inds.first].second.push(subrun, trigger_offset);
		}

		bool contains(const int& run, const int& subrun) const
		{
			const indexed_info<T> *subruns = ((indexed_info<indexed_info<T> >*) this)->info(run);
			if (NULL == subruns)
				return false;
			return subruns->contains(subrun);
		}

		T operator() (const int& run, const int& subrun) const //returns std::numeric_limits<T>::max() if not found
		{
			const indexed_info<T> *subruns = ((indexed_info<indexed_info<T> >*) this)->info(run);
			if (NULL == subruns)
				return std::numeric_limits<T>::max();
			const T *val = subruns->info(subrun);
			if (NULL == val)
				return std::numeric_limits<T>::max();
			return *val;
		}

		T* info(const int& run, const int& subrun)
		{
			indexed_info<T> *subruns = ((indexed_info<indexed_info<T> >*) this)->info(run);
			if (NULL == subruns)
				return NULL;
			return subruns->info(subrun);
		}
		const T* info(const int& run, const int& subrun) const
		{
			const indexed_info<T> *subruns = ((indexed_info<indexed_info<T> >*) this)->info(run);
			if (NULL == subruns)
				return NULL;
			return subruns->info(subrun);
		}
	};

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

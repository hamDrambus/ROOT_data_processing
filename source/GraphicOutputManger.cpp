#include "GraphicOutputManager.h"

Drawing::Drawing(std::string name, ParameterPile::DrawEngine de, int id_index):
_name(name), _id_index(id_index), _de(de), _directory(OUTPUT_DIR+"gnuplot/")
{
	for (int s = 0; s < _name.size(); s++)
		if (_name[s] == '\\' || _name[s] == '/')
			_name[s] = '.';
	_script_fname = _name + ".sc";
	_script_lines.push_back("clear");
	_script_lines.push_back("set multiplot");
	_script_lines.push_back("set ytics nomirror");
	_script_lines.push_back("set y2tics");
	_script_lines.push_back("set key top right");
	_script_lines.push_back("###"); //pad-specific commands are between ### lines
	_script_lines.push_back("###");
	_script_lines.push_back("unset multiplot");
#if !defined(__WIN32__)
	_script_lines.push_back("pause -1");
#endif
}

int Drawing::get_index_of_pad_marker(int pad_index)//before which line enter the script;
{
	if (pad_index < 0)
		return -1;

	int marker_n = 0;
	int off = 0;
	int last_pos = 0;
	for (auto i = _script_lines.begin(); (i != _script_lines.end()) && (marker_n < (pad_index + 2)); i++, off++)
		if ((*i) == std::string("###")){
			marker_n++;
			last_pos = off;
		}
	last_pos++;
	off--;//before which line enter the script;
	while (marker_n < (pad_index + 2)){
		_script_lines.insert(_script_lines.begin() + last_pos, "###");
		marker_n++;
		last_pos++;
		off = last_pos - 1;
	}
	return off;
}

std::string Drawing::process_title (std::string in)
{
	std::string out;
	char prev_ch = 'a';
	for (std::size_t i=0, i_end_=in.size();i!=i_end_;++i) {
		if ((in[i]=='_')&&(prev_ch!='\\')) {
			out.push_back('\\');
			out.push_back(in[i]);
			prev_ch = in[i];
		} else {
			out.push_back(in[i]);
			prev_ch = in[i];
		}
	}
	return out;
}

void Drawing::AddToDraw(DVECTOR &xs, DVECTOR &ys, std::string title, std::string extra_txt, int pad_index)
{
	int off = get_index_of_pad_marker(pad_index);
	if (off < 0)
		return;
	title = process_title(title);
	_data_fnames.push_back(_directory+_name +"_"+ std::to_string(_data_fnames.size()));
	std::string line = "plot '"+ParameterPile::this_path + _data_fnames.back() +"' u 1:2 title '"+title+"' " + extra_txt;
	_script_lines.insert(_script_lines.begin() + off, line);

	line = _data_fnames.back();
	std::ofstream f_data;
	open_output_file(line, f_data);
	for (auto i = xs.begin(), j = ys.begin(); (i != xs.end()) && (j != ys.end()); i++, j++)
		f_data<< *i << '\t' << *j << std::endl;
	f_data.close();
}

void Drawing::AddToDraw(DVECTOR &xs, DVECTOR &ys, DVECTOR &ys_err, std::string title, std::string extra_txt, int pad_index)
{
	int off = get_index_of_pad_marker(pad_index);
	if (off < 0)
		return;
	title = process_title(title);
	_data_fnames.push_back(_directory+_name +"_"+ std::to_string(_data_fnames.size()));
	std::string line = "plot '"+ParameterPile::this_path + _data_fnames.back() +"' u 1:2:3 with errorbars title '"+title+"' " + extra_txt;
	_script_lines.insert(_script_lines.begin() + off, line);

	line = _data_fnames.back();
	std::ofstream f_data;
	open_output_file(line, f_data);
	for (auto i = xs.begin(), j = ys.begin(), e = ys_err.begin(), i_end_ = xs.end(), j_end_=ys.end(), e_end_ = ys_err.end();
			(i != i_end_) && (j != j_end_) && (e != e_end_); ++i, ++j, ++e)
		f_data<< *i << '\t' << *j << '\t' << *e <<std::endl;
	f_data.close();
}

void Drawing::AddToDraw_baseline(double base_line, std::string title, std::string extra_txt, int pad_index)
{
	std::string definition_lines;
	AddToDraw(definition_lines, std::to_string(base_line), title, extra_txt, pad_index);
}

void Drawing::AddToDraw_vertical(double x_pos, double from_y, double to_y, std::string extra_txt, int pad_index)
{
	int off = get_index_of_pad_marker(pad_index);
	if (off < 0)
		return;
	std::vector<std::string> new_lines;
	new_lines.push_back("set arrow from "+std::to_string(x_pos)+","+std::to_string(from_y)+" to "+std::to_string(x_pos)+","+std::to_string(to_y)+" nohead "
		+extra_txt);
	new_lines.push_back("show arrow");
	_script_lines.insert(_script_lines.begin() + off, new_lines.begin(), new_lines.end());
}

void Drawing::AddToDraw(std::string definition_lines, std::string f_name, std::string title, std::string extra_txt, int pad_index)//any gnuplot function
{
	int off = get_index_of_pad_marker(pad_index);
	if (off < 0)
		return;
	STD_CONT<std::string> new_lines;
	title = process_title(title);
	new_lines.push_back(definition_lines);
	new_lines.push_back("plot "+f_name+" title '"+title+"' "+extra_txt);
	_script_lines.insert(_script_lines.begin() + off, new_lines.begin(), new_lines.end());
}

void Drawing::DrawData(DVECTOR &xs, DVECTOR &ys, std::string title, std::string extra_txt)//draws only this vector
{
	if (xs.size() != ys.size()){
		std::cout << "Drawing::DrawData::input data size mismatch" << std::endl;
		return;
	}
	if (_de == ParameterPile::DrawEngine::ROOT){
		double *xxxs = new double[xs.size()];
		double *yyys = new double[ys.size()];
		for (int h = 0; h < xs.size(); h++){
			xxxs[h] = xs.at(h);
			yyys[h] = ys.at(h);
		}
		TGraph* gr = new TGraph(xs.size(), xxxs, yyys);
		TCanvas* can = new TCanvas(_name.c_str(), title.c_str());
		can->SetWindowPosition(100, 150);
		can->Update();
		gr->Draw();
		delete[] xxxs;
		delete[] yyys;
	} else {
		std::string mod_name = _name;
		for (int s = 0; s < mod_name.size(); s++)
			if (mod_name[s] == '\\' || mod_name[s] == '/')
				mod_name[s] = '.';
		std::ofstream file;
		open_output_file(_directory + mod_name, file);
		std::cout << "file '" << _directory + mod_name << "'.is_open() " << file.is_open() << std::endl;
		if (!file.is_open()){
			std::cout<<"Could not open a file!";
#if defined (__WIN32__)
			std::cout << GetLastError() << std::endl;
#endif //__WIn32__
		}
		for (int h = 0; h < xs.size(); h++)
			file << xs[h] << '\t' << ys[h] << std::endl;
		file.close();
		open_output_file(_directory +_script_fname, file);
		title = process_title(title);
		file << "plot '" << ParameterPile::this_path + _directory + mod_name << "' u 1:2 title '" << title<<"' "<<extra_txt<< std::endl;
		file << "pause -1";
		file.close();
		INVOKE_GNUPLOT(ParameterPile::this_path + _directory + _script_fname);
	}
}

void Drawing::DrawData(void)
{
	//int N_marks =0;
	int N_pads = 0;
	bool is_empty = false;
	//this is set to false in order to account properly for the first "###"
	for (auto i = _script_lines.begin(); i != _script_lines.end(); i++){
		if ((*i) == "###"){
			if (!is_empty)
				N_pads++;
			is_empty = true;
		} else 
			is_empty = (*i).empty() ? is_empty : false;
	}
	N_pads--; //must still be >=1, so no check here
	if (1 == N_pads){//no need in multiplot
		_script_lines[1].insert(0,"#");
	}

#if defined(__WIN32__)
	_script_lines.insert(_script_lines.begin() + 1, "set terminal wxt size " + std::to_string(ParameterPile::gnuplot_width) + ","
		+ std::to_string(std::min(ParameterPile::gnuplot_max_size , N_pads*ParameterPile::gnuplot_pad_size)));
#else
	_script_lines.insert(_script_lines.begin() + 1, "set terminal qt size " + std::to_string(ParameterPile::gnuplot_width) + ","
			+ std::to_string(std::min(ParameterPile::gnuplot_max_size , N_pads*ParameterPile::gnuplot_pad_size)));
#endif
	int set_pads = 0;
	while (set_pads < N_pads){ //because inserting the lines invalidates iterators
		int pads_off = set_pads + 1;//+1 in order to account for the first "###"
		for (auto i = _script_lines.begin(); i != _script_lines.end(); i++)
			if ((*i) == std::string("###")){
				is_empty = true;
				bool has_plot = false;
				for (auto j = (i + 1); (j != _script_lines.end()) && ((*j) != "###"); ++j) {//scan the pad
					//There should be only one 'plot' in one pad and the rest are 'replot'. Otherwise the same color is used and title is simply overriden 
					if (0 == (*j).find("plot")){
						if (has_plot)
							(*j).insert(0, "set yrange restore\nre");
						else {
							(*j).insert(0, "set yrange [] writeback\n"); //so my y scale in each pad is defined by the first added function
							has_plot = true;
						}
					}
					is_empty = (*j).empty() ? is_empty : false;
				}
				if (!is_empty){
					pads_off--;
					if (0 == pads_off){
						STD_CONT<std::string> pad_lines;
						pad_lines.push_back("unset arrow");
						pad_lines.push_back("set origin 0,"+std::to_string(set_pads*(1.0/N_pads)));
						pad_lines.push_back("set size 1," + std::to_string(1.0 / N_pads));
						_script_lines.insert(i + 1, pad_lines.begin(), pad_lines.end());
						set_pads++;
						break;
					}
				}
			}
	}
	std::ofstream file;
	open_output_file(_directory + _script_fname, file);
	for (auto i = _script_lines.begin(); i != _script_lines.end(); i++)
		file << *i << std::endl;
	file.close();
	INVOKE_GNUPLOT(ParameterPile::this_path + _directory + _script_fname);
}

void Drawing::Clear(void)
{

}

void Drawing::SetDirectory(std::string path)
{
	_directory = path;
}

std::string Drawing::get_name(void) const
{	return _name;}
int Drawing::get_id(void) const
{	return _id_index;}

GraphicOutputManager::GraphicOutputManager(void)
{}

Drawing* GraphicOutputManager::GetDrawing(int index) //if does not exist, doest not create it
{
	for (auto i = _graphs.begin(); i != _graphs.end(); i++)
		if ((*i).get_id() == index)
			return &(*i);
	return NULL;
}

Drawing* GraphicOutputManager::GetDrawing(std::string name, int index, ParameterPile::DrawEngine de)//if does not exist, creates it
{
	for (auto i = _graphs.begin(); i != _graphs.end(); i++)
		if ((*i).get_name() == name || (*i).get_id() == index)
			return &(*i);
	return CreateDrawing(name, index, de);
}

Drawing* GraphicOutputManager::CreateDrawing(std::string name, int index, ParameterPile::DrawEngine de)
{
	for (auto i = _graphs.begin(); i != _graphs.end(); i++)
		if ((*i).get_name() == name || (*i).get_id() == index)
			return NULL;
	_graphs.push_back(Drawing(name, de, index));
	return &_graphs.back();
}

void GraphicOutputManager::Draw(void)
{
	for (auto i = _graphs.begin(); i != _graphs.end(); i++)
		(*i).DrawData();
}

void GraphicOutputManager::Clear(void)
{
#ifdef _HOTFIX_CLEAR_MEMORY
	STD_CONT<Drawing>().swap(_graphs);
#else
	_graphs.clear();
#endif
}

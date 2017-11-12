#ifndef GRAPHIC_OUTPUT_MANAGER_H
#define GRAPHIC_OUTPUT_MANAGER_H

#include "GlobalParameters.h"

class Drawing
{
protected:
	std::string _name;
	int _id_index;
	ParameterPile::DrawEngine _de;
	STD_CONT<std::string> _data_fnames;
	std::string _script_fname;
	STD_CONT<std::string> _script_lines;

	int get_index_of_pad_marker(int pad);//before which line enter the script;

public:
	//Drawing(void);
	Drawing(std::string name, ParameterPile::DrawEngine de, int id_index);//can't change DrawEngine - too much of a trouble
	void AddToDraw(DVECTOR &xs, DVECTOR &ys, std::string title = "", std::string extra_txt = "", int pad_index = 0);
	void AddToDraw_baseline(double base_line, std::string title = "baseline", std::string extra_txt = "", int pad_index = 0);//May add more functions to draw, e.g. gauss
	void AddToDraw_vertical(double x_pos, double from_y, double to_y, std::string extra_txt="", int pad_index = 0);
	void AddToDraw(std::string definition_lines, std::string f_name, std::string title, std::string extra_txt = "", int pad_index = 0);//any gnuplot function
	
	void DrawData(DVECTOR &xs, DVECTOR &ys, std::string title = "", std::string extra_txt = "");//draws only this vector
	void DrawData(void);
	void Clear(void);

	std::string get_name(void) const;
	int get_id(void) const;
};

class GraphicOutputManager
{
public:
	STD_CONT<Drawing> _graphs;
	GraphicOutputManager(void);
	Drawing* GetDrawing(int index); //if does not exist, doest not create it
	Drawing* GetDrawing(std::string name, int index, ParameterPile::DrawEngine de);//if does not exist, creates it
	Drawing* CreateDrawing(std::string name, int index, ParameterPile::DrawEngine de);
	void Draw (void);
	void Clear (void);
};

#endif
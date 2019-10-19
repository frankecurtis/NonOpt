// Copyright (C) 2019 Minhan Li
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#ifndef __SETDIM_HPP__
#define __SETDIM_HPP__


class setDim
{

 public:
setDim(){
	 dim_=50;
 };

 ~setDim(){};

 inline int getDim() const {return dim_;};

 private:
 int dim_;
};

#endif

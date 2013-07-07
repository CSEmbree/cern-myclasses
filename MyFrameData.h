#ifndef __MY_FRAME_DATA_H
#define __MY_FRAME_DATA_H

#include "root.h"
#include "MyFrame.h"
#include "MyData.h"

class MyFrameData : MyFrame 
//class MyFrameData : 
{
 private:

  bool debug;
  MyFrame *myframe ;

 public:

  MyFrameData(Int_t p1, Int_t p2, MyData *gdata);
  MyFrame * GetMyFrame(){ return myframe;}  
  //void Print();
  void SaveFile(TString file_name);

};
#endif


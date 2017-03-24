
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _AMReX_EBDataVarMacros_H_
#define _AMReX_EBDataVarMacros_H_

#if (BL_SPACEDIM == 2)

#define V_VOLFRAC       0
#define V_BNDAREA       1
#define V_VOLCENTROIDX  2
#define V_VOLCENTROIDY  3
#define V_BNDCENTROIDX  4
#define V_BNDCENTROIDY  5
#define V_NORMALX       6
#define V_NORMALY       7
#define V_VOLNUMBER     8
                        
#define F_AREAFRAC      0
#define F_FACECENTROIDX 1
#define F_FACECENTROIDY 2
#define F_FACENUMBER    3

#elif (BL_SPACEDIM == 3)

#define V_VOLFRAC      0
#define V_BNDAREA      1
#define V_VOLCENTROIDX 2
#define V_VOLCENTROIDY 3
#define V_VOLCENTROIDZ 4
#define V_BNDCENTROIDX 5
#define V_BNDCENTROIDY 6
#define V_BNDCENTROIDZ 7
#define V_NORMALX      8
#define V_NORMALY      9
#define V_NORMALZ      10
#define V_VOLNUMBER    11

#define F_AREAFRAC      0
#define F_FACECENTROIDX 1
#define F_FACECENTROIDY 2
#define F_FACECENTROIDZ 3
#define F_FACENUMBER    4

#endif


#endif

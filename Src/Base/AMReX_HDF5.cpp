#include <AMReX_HDF5.H>
#include <iostream>

#ifdef AMREX_USE_HDF5

namespace amrex
{

hid_t makeH5Box() {
  hid_t box_id = H5Tcreate(H5T_COMPOUND, sizeof(box_h5_t));
#if AMREX_SPACEDIM >= 1
  H5Tinsert(box_id, "lo_i", HOFFSET(box_h5_t, lo_i), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(box_id, "lo_j", HOFFSET(box_h5_t, lo_j), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(box_id, "lo_k", HOFFSET(box_h5_t, lo_k), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 1
  H5Tinsert(box_id, "hi_i", HOFFSET(box_h5_t, hi_i), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(box_id, "hi_j", HOFFSET(box_h5_t, hi_j), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(box_id, "hi_k", HOFFSET(box_h5_t, hi_k), H5T_NATIVE_INT);
#endif

  return box_id;
}

box_h5_t writeH5Box(const Box &b) {
  box_h5_t box;
#if AMREX_SPACEDIM >= 1
  box.lo_i = b.smallEnd(0);
  box.hi_i = b.bigEnd(0);
#endif
#if AMREX_SPACEDIM >= 2
  box.lo_j = b.smallEnd(1);
  box.hi_j = b.bigEnd(1);
#endif
#if AMREX_SPACEDIM >= 3
  box.lo_k = b.smallEnd(2);
  box.hi_k = b.bigEnd(2);
#endif
  return box;
}

void writeH5Box(const Box &b, box_h5_t &box) {
#if AMREX_SPACEDIM >= 1
  box.lo_i = b.smallEnd(0);
  box.hi_i = b.bigEnd(0);
#endif
#if AMREX_SPACEDIM >= 2
  box.lo_j = b.smallEnd(1);
  box.hi_j = b.bigEnd(1);
#endif
#if AMREX_SPACEDIM >= 3
  box.lo_k = b.smallEnd(2);
  box.hi_k = b.bigEnd(2);
#endif
  return;
}

Box readH5Box(box_h5_t &box) {
  IntVect lo(AMREX_D_DECL(box.lo_i, box.lo_j, box.lo_k));
  IntVect hi(AMREX_D_DECL(box.hi_i, box.hi_j, box.hi_k));
  Box b(lo, hi);
  return b;
}

void writeBoxOnHDF5(const Box& box, H5& h5, const std::string name)
{
  hid_t box_id = makeH5Box();
  box_h5_t h5box = writeH5Box(box);
  h5.writeAttribute(name, h5box, box_id);
  H5Tclose(box_id);
  return;
}

Box readBoxFromHDF5(H5& h5, const std::string name)
{
  box_h5_t box;
  h5.readAttribute(name, box);
  Box out = readH5Box(box);
  return out;
}


hid_t makeH5RealBox() {
  hid_t box_id = H5Tcreate(H5T_COMPOUND, sizeof(rbox_h5_t));
#if AMREX_SPACEDIM >= 1
  H5Tinsert(box_id, "lo_x", HOFFSET(rbox_h5_t, lo_x), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(box_id, "lo_y", HOFFSET(rbox_h5_t, lo_y), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(box_id, "lo_z", HOFFSET(rbox_h5_t, lo_z), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 1
  H5Tinsert(box_id, "hi_x", HOFFSET(rbox_h5_t, hi_x), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(box_id, "hi_y", HOFFSET(rbox_h5_t, hi_y), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(box_id, "hi_z", HOFFSET(rbox_h5_t, hi_z), H5T_NATIVE_DOUBLE);
#endif

  return box_id;
}

rbox_h5_t writeH5RealBox(const RealBox &b) {
  rbox_h5_t box;
#if AMREX_SPACEDIM >= 1
  box.lo_x = b.lo(0);
  box.hi_x = b.hi(0);
#endif
#if AMREX_SPACEDIM >= 2
  box.lo_y = b.lo(1);
  box.hi_y = b.hi(1);
#endif
#if AMREX_SPACEDIM >= 3
  box.lo_z = b.lo(2);
  box.hi_z = b.hi(2);
#endif
  return box;
}

RealBox readH5RealBox(rbox_h5_t &box) {
  std::array<Real,AMREX_SPACEDIM> lo = {AMREX_D_DECL(box.lo_x, box.lo_y, box.lo_z)};
  std::array<Real, AMREX_SPACEDIM> hi = {AMREX_D_DECL(box.hi_x, box.hi_y, box.hi_z)};
  RealBox b(lo, hi);
  return b;
}

void writeRealBoxOnHDF5(const RealBox& box, H5& h5, const std::string name)
{
  hid_t rbox_id = makeH5RealBox();
  rbox_h5_t rbox = writeH5RealBox(box);
  h5.writeAttribute(name, rbox, rbox_id);
  H5Tclose(rbox_id);
  return;
}

RealBox readRealBoxFromHDF5(H5& h5, const std::string name)
{
  rbox_h5_t rbox;
  h5.readAttribute(name, rbox);
  RealBox out = readH5RealBox(rbox);
  return out;
}

hid_t makeH5IntVec() {
  hid_t intvect_id = H5Tcreate(H5T_COMPOUND, sizeof(int_h5_t));
#if AMREX_SPACEDIM >= 1
  H5Tinsert(intvect_id, "intvecti", HOFFSET(int_h5_t, i), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(intvect_id, "intvectj", HOFFSET(int_h5_t, j), H5T_NATIVE_INT);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(intvect_id, "intvectk", HOFFSET(int_h5_t, k), H5T_NATIVE_INT);
#endif
  return intvect_id;
}

int_h5_t writeH5IntVec(const int* in) {
  int_h5_t i;
#if AMREX_SPACEDIM >= 1
  i.i = in[0];
#endif
#if AMREX_SPACEDIM >= 2
  i.j = in[1];
#endif
#if AMREX_SPACEDIM >= 3
  i.k = in[2];
#endif
  return i;
}

void readH5IntVec(int_h5_t &in, int *out) {
#if AMREX_SPACEDIM >= 1
  out[0] = in.i;
#endif
#if AMREX_SPACEDIM >= 2
  out[1] = in.j;
#endif
#if AMREX_SPACEDIM >= 3
  out[2] = in.k;
#endif
}

hid_t makeH5RealVec() {
hid_t realvect_id = H5Tcreate(H5T_COMPOUND, sizeof(real_h5_t));
#if AMREX_SPACEDIM >= 1
  H5Tinsert(realvect_id, "x", HOFFSET(real_h5_t, x), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 2
  H5Tinsert(realvect_id, "y", HOFFSET(real_h5_t, y), H5T_NATIVE_DOUBLE);
#endif
#if AMREX_SPACEDIM >= 3
  H5Tinsert(realvect_id, "z", HOFFSET(real_h5_t, z), H5T_NATIVE_DOUBLE);
#endif
  return realvect_id;
}

real_h5_t writeH5RealVec(const Real *in) {
  real_h5_t vec;
#if AMREX_SPACEDIM >= 1
  vec.x = in[0];
#endif
#if AMREX_SPACEDIM >= 2
  vec.y = in[1];
#endif
#if AMREX_SPACEDIM >= 3
  vec.z = in[2];
#endif
  return vec;
}

void readH5RealVec(real_h5_t &in, double *out) {
#if AMREX_SPACEDIM >= 1
  out[0] = in.x;
#endif
#if AMREX_SPACEDIM >= 2
  out[1] = in.y;
#endif
#if AMREX_SPACEDIM >= 3
  out[2] = in.z;
#endif
}

//==================================================================================
// H5
//==================================================================================


H5::H5() {}

H5::H5(std::string name) { createFile(name); }

H5::H5(std::string name, MPI_Comm comm) { createFile(name, comm); }

H5::H5(hid_t h5) { m_obj = h5; }

H5::~H5() {}

void H5::createFile(const std::string name, MPI_Comm comm) {
  // assume MPI already initialized
  hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_access, comm, MPI_INFO_NULL);

  m_obj = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_access);
  m_name = name;
  H5Pclose(file_access);

  track(name);

  return;
}

void H5::openFile(const std::string name, MPI_Comm comm) {


  hid_t file_access = H5Pcreate (H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_access,  comm, MPI_INFO_NULL);


  m_obj = H5Fopen(name.c_str(), H5F_ACC_RDWR, file_access);
  if (m_obj < 0) {
    amrex::Abort("Unbale to open " + name);
  }

  m_name = name;

  H5Pclose(file_access);

  track(name);

  return;
}

void H5::closeFile() {
  discard(m_name);
#ifdef AMREX_DEBUG
  if (!tracker.empty()) {
    std::cout << "H5 tracker outstanding items";
    for (auto& name : tracker) {
      std::cout << name << "\n";
    }
  }
#endif

  H5Fclose(m_obj);

}

bool H5::groupExists(const std::string name) {
  const char* msg="an error has occured\n";
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, (void*)msg);
  status = H5Gget_objinfo (m_obj, name.c_str(), 0, NULL);

  if (status == 0) {
    return true;
  } else {
    return false;
  }
}

H5 H5::createGroup(const std::string name) {
  H5 out;

  if (groupExists(name)) {
    out.m_obj = H5Gopen2(m_obj, name.c_str(), H5P_DEFAULT);
  } else {
    out.m_obj = H5Gcreate(m_obj, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  if (out.m_obj < 0) {
    amrex::Abort(" Problem creating group " + name);
  }
  out.m_name = name;
  track(name);
  return out;
}

H5 H5::openGroup(const std::string name) {
  H5 out;
  out.m_obj = H5Gopen2(m_obj, name.c_str(), H5P_DEFAULT);
  if (out.m_obj < 0) {
    amrex::Abort(" Problem opening group " + name);
  }
  out.m_name = name;
  track(name);
  return out;
}

void H5::closeGroup() {
  H5Gclose(m_obj);
  discard(m_name);
}

H5 H5::openDataset(const std::string name) {
  H5 out;
  out.m_obj = H5Dopen2(m_obj, name.c_str(), H5P_DEFAULT);
  out.m_name = name;
  track(name);
  return out;
}

void H5::closeDataset() {
  H5Dclose(m_obj);
  discard(m_name);
}

herr_t H5::writeAttribute(std::map<std::string, int>& m_int,
                         std::map<std::string, double>& m_real,
                         std::map<std::string, std::string>& m_string) {
  H5E_auto_t efunc;
  void* edata;
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
  herr_t ret;

#define INSERT2(Ttype, mapName, H5Ttype)                                      \
  for (std::map<std::string, Ttype>::const_iterator p = mapName.begin();      \
       p != mapName.end(); ++p) {                                             \
    hid_t aid = H5Screate(H5S_SCALAR);                                        \
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);                                    \
    hid_t attr = H5Acreate2(m_obj, p->first.c_str(), H5Ttype, aid, H5P_DEFAULT, \
                            H5P_DEFAULT);                                     \
    if (attr < 0) {                                                           \
      H5Adelete(m_obj, p->first.c_str());                                       \
      attr = H5Acreate2(m_obj, p->first.c_str(), H5Ttype, aid, H5P_DEFAULT,     \
                        H5P_DEFAULT);                                         \
      if (attr < 0) {                                                         \
        amrex::Abort(" Problem writing attribute " + p->first);       \
      }                                                                       \
    }                                                                         \
    H5Eset_auto2(H5E_DEFAULT, efunc, edata);                                  \
    Ttype tmp = p->second;                                                    \
    ret = H5Awrite(attr, H5Ttype, &tmp);                                      \
    if (ret < 0) return ret;                                                  \
    H5Sclose(aid);                                                            \
    H5Aclose(attr);                                                           \
  }
  INSERT2(double, m_real, H5T_NATIVE_DOUBLE);
  INSERT2(int, m_int, H5T_NATIVE_INT);

  // string is different, of course
  for (std::map<std::string, std::string>::const_iterator p = m_string.begin();
       p != m_string.end(); ++p) {
    hid_t s_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(s_type, p->second.length());  // extra requirement for strings
    hid_t aid = H5Screate(H5S_SCALAR);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    hid_t attr = H5Acreate2(m_obj, p->first.c_str(), s_type, aid, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (attr < 0) {
      H5Adelete(m_obj, p->first.c_str());
      attr = H5Acreate2(m_obj, p->first.c_str(), s_type, aid, H5P_DEFAULT,
                        H5P_DEFAULT);
      if (attr < 0) {
        amrex::Abort(" Problem writing attribute " + p->first);
      }
    }
    H5Eset_auto2(H5E_DEFAULT, efunc, edata);
    char* tmp = (char*)p->second.c_str();
    ret = H5Awrite(attr, s_type, tmp);
    if (ret < 0) return ret;
    H5Sclose(aid);
    H5Aclose(attr);
    H5Tclose(s_type);
  }

  return 0;
}

void H5::writeString(const std::string name, const std::string& data) {
  // write a single string to obj

  hid_t type, space, dset;

  type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, data.size());

  hsize_t dims[1] = {1};
  space = H5Screate_simple(1, dims, NULL);

  dset = H5Dcreate(m_obj, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT,
                   H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(type);

  return;
}

void H5::writeString(const std::string name,
                    const std::vector<std::string>& data) {
  // write a list of strings
  // all strings must be the same length

  // create a 1d buffer of all the characters

  std::vector<char> buffer;

  for (std::string S : data) {
    for (size_t i = 0; i < S.size(); ++i) {
      buffer.push_back(S[i]);
    }
  }

  hsize_t dims[1] = {data.size()};
  hid_t type, space, dset;

  type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, data[0].size());

  space = H5Screate_simple(1, dims, NULL);

  dset = H5Dcreate(m_obj, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT,
                   H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(type);

  return;
}

std::list<std::string> H5::tracker = {};
void H5::track(const std::string name)
{
#ifdef AMREX_DEBUG
  tracker.push_back(name);
#endif
}

void H5::discard(const std::string name)
{
#ifdef AMREX_DEBUG
  tracker.remove(name);
#endif
}


}

#endif

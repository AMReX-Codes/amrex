
#include <iostream.h>
#include <malloc.h>
#include "amr_defs.H"

void malloc_info()
{
#ifdef BL_Linux
  struct mstats mi = mstats();
  cout << "Total Heap Size: " << mi.bytes_total << endl;
  cout << "Allocated: " << mi.chunks_used << " " << mi.bytes_used << endl;
  cout << "Free:      " << mi.chunks_free << " " << mi.bytes_free << endl;
#else
  struct mallinfo mi = mallinfo();

  cout << "Malloc arena: " << mi.arena << endl;
  cout << "Ordinary blocks: " << mi.ordblks << " "
    << mi.uordblks << " " << mi.fordblks << endl;
#  ifndef SYSV
  cout << "Allocated: " << mi.allocated << " " << mi.uordbytes << endl;
#  endif

  if (mi.smblks || mi.usmblks || mi.fsmblks || mi.hblks || mi.hblkhd) {
    cout << "Small blocks: "
      << mi.smblks << " " << mi.usmblks << " " << mi.fsmblks << " "
      << mi.hblks << " " << mi.hblkhd << endl;
  }
#endif
}

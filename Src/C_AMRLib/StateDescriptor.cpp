//BL_COPYRIGHT_NOTICE

//
// $Id: StateDescriptor.cpp,v 1.11 1999-05-10 18:54:08 car Exp $
//

#include <StateDescriptor.H>
#include <Interpolater.H>
#include <BCRec.H>

DescriptorList::DescriptorList ()
    :
    desc(PArrayManage)
{}

DescriptorList::~DescriptorList () {}

void
DescriptorList::addDescriptor (int                         indx,
                               IndexType                   typ,
                               StateDescriptor::TimeCenter ttyp,
                               int                         nextra,
                               int                         num_comp, 
                               Interpolater*               interp)
{
    if (indx >= desc.length())
        desc.resize(indx+1);
    desc.set(indx, new StateDescriptor(typ,ttyp,indx,nextra,num_comp,interp));
}  

StateDescriptor::StateDescriptor ()
    :
    id(-1),
    ncomp(0),
    ngrow(0),
    mapper(0),
    bc_func(PArrayManage),
    t_type(Point)
{}

StateDescriptor::StateDescriptor (IndexType                   btyp,
                                  StateDescriptor::TimeCenter ttyp,
                                  int                         ident,
                                  int                         nextra, 
                                  int                         num_comp,
                                  Interpolater*               interp)
    :
    type(btyp),
    t_type(ttyp),
    id(ident),
    ngrow(nextra),
    ncomp(num_comp),
    mapper(interp),
    bc_func(PArrayManage)
{
    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

StateDescriptor::~StateDescriptor ()
{
    mapper = 0;
}

void
StateDescriptor::define (IndexType                   btyp,
                         StateDescriptor::TimeCenter ttyp,
                         int                         ident,
                         int                         nextra,
                         int                         num_comp,
                         Interpolater*               interp) 
{
    type   = btyp;
    t_type = ttyp;
    id     = ident;
    ngrow  = nextra;
    ncomp  = num_comp;
    mapper = interp;

    BL_ASSERT (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

void
StateDescriptor::setComponent (int                               comp,
                               const aString&                    nm,
                               const BCRec&                      bcr,
                               const StateDescriptor::BndryFunc& func,
                               Interpolater*                     interp, 
                               int                               max_map_start_comp_,
                               int                               min_map_end_comp_)
{
    BL_ASSERT(comp >= 0 && comp < ncomp && names[comp].isNull());
    names[comp]       = nm;
    bc_func.set(comp,func.clone());
    bc[comp]          = bcr;
    mapper_comp[comp] = interp;

    if (max_map_start_comp_>=0 && min_map_end_comp_>=0)
    {
        BL_ASSERT(comp >= max_map_start_comp_ &&
               comp <= min_map_end_comp_   &&
               min_map_end_comp_ < ncomp);
        max_map_start_comp[comp] = max_map_start_comp_;
        min_map_end_comp[comp]   = min_map_end_comp_;
    }
    else
    {
        max_map_start_comp[comp] = comp;
        min_map_end_comp[comp]   = comp;
    }
}

void
StateDescriptor::dumpNames (ostream& os,
                            int      start_comp,
                            int      num_comp) const
{
    BL_ASSERT(start_comp >= 0 && start_comp+num_comp <= ncomp);

    for (int k = 0; k < num_comp; k++)
    {
        os << names[start_comp+k] << ' ';
    }
}

void
StateDescriptor::setUpMaps (int&                use_default_map,
                            const Interpolater* default_map,
                            int                 start_comp,
                            int                 num_comp,
                            Interpolater**&     maps, 
                            int&                nmaps,
                            int*&               map_start_comp,
                            int*&               map_num_comp, 
                            int*&               max_start_comp,
                            int*&               min_end_comp) const
{
    BL_ASSERT(start_comp>=0 && start_comp+num_comp-1 < ncomp && num_comp>0);

    maps           = 0;
    map_start_comp = 0;
    map_num_comp   = 0;
    max_start_comp = 0;
    min_end_comp   = 0;
    //
    // First, count number of interpolaters needed and allocate.
    //
    Interpolater* map = mapper_comp[start_comp];
    if (!map) map = (Interpolater*) default_map;
    nmaps = 1; 
    int icomp = start_comp+1;

    use_default_map = 1;
    while (icomp < start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];
        if (!mapper_icomp)
        {
            mapper_icomp = (Interpolater *) default_map;
        }
        else
        {
            use_default_map = 0;
        }
        if (map != mapper_icomp)
        {
            map = mapper_icomp;
            nmaps++;
        }
        icomp++;
    }

    if (use_default_map) return;

    maps           = new Interpolater*[nmaps];
    map_start_comp = new int[nmaps];
    map_num_comp   = new int[nmaps];
    min_end_comp   = new int[nmaps];
    max_start_comp = new int[nmaps];
    //
    // Now fill the slots.
    //
    int imap             = 0;
    if (mapper_comp[start_comp])
        maps[imap]         = mapper_comp[start_comp];
    else 
        maps[imap]         = (Interpolater *) default_map;

    icomp                = start_comp+1;
    map_start_comp[imap] = start_comp;
    map_num_comp[imap]   = 1;
    
    min_end_comp[imap]   = min_map_end_comp[start_comp];
    max_start_comp[imap] = max_map_start_comp[start_comp];

    while (icomp<start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];
        if (!mapper_icomp) mapper_icomp = (Interpolater *) default_map;

        if (maps[imap] != mapper_icomp)
        {
            imap++;

            BL_ASSERT (imap < nmaps);

            maps[imap]           = mapper_icomp;
            map_start_comp[imap] = icomp;
            map_num_comp[imap]   = 1;
            min_end_comp[imap]   = min_map_end_comp[icomp];
            max_start_comp[imap] = max_map_start_comp[icomp];

        }
        else
        {
            map_num_comp[imap]++;
            min_end_comp[imap]   = Max(min_end_comp[imap],min_map_end_comp[icomp]);
            max_start_comp[imap] = Min(max_start_comp[imap],max_map_start_comp[icomp]);
        }
        icomp++;
    }
}

void
StateDescriptor::cleanUpMaps (Interpolater**& maps, 
                              int*&           map_start_comp,
                              int*&           map_num_comp,
                              int*&           max_start_comp,
                              int*&           min_end_comp) const
{
    delete maps;
    delete map_start_comp;
    delete map_num_comp;
    delete max_start_comp;
    delete min_end_comp;
}

// $Id:$

#include <stdio.h>
#include "CMS2TreePlayer.h"
#include "CMS2Selector.h"

bool swap_treeplayer ();
static bool treeplayer_swapped_ = swap_treeplayer();

bool swap_treeplayer ()
{
     TVirtualTreePlayer::SetPlayer("CMS2TreePlayer");
     return true;
}

CMS2TreePlayer::CMS2TreePlayer ()
     : TTreePlayer() 
{ 
     printf("Oooh, I'm a CMS2TreePlayer!\n"); 
     fSelector = new CMS2Selector;
     fSelector->SetInputList(fInput);
}

ClassImp(CMS2TreePlayer)

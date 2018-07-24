/*
 * ntuple_global.h
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_

class ntuple_global: public ntuple_content{
public:
	ntuple_global():ntuple_content(){clear();}

private:
	float numPileUp;


};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_ */

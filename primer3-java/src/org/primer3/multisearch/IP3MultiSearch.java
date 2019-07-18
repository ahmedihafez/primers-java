package org.primer3.multisearch;

import org.primer3.libprimer3.DPAlArgHolder;
import org.primer3.libprimer3.THAlArgHolder;

public interface IP3MultiSearch {

	DPAlArgHolder getDPAL_ArgToUse();

	THAlArgHolder get_thal_primers_arg_to_use();

}

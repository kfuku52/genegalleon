"""Compatibility CLI for NWKIT calibration unit conversion."""

from decimal import Decimal as Decimal

from nwkit.calibration_scale import CALIBRATION_RE as CALIBRATION_RE
from nwkit.calibration_scale import NUMERIC_RE as NUMERIC_RE
from nwkit.calibration_scale import (
    build_parser as build_parser,
)
from nwkit.calibration_scale import (
    calibration_age_field_count as calibration_age_field_count,
)
from nwkit.calibration_scale import (
    calibration_ages as calibration_ages,
)
from nwkit.calibration_scale import (
    choose_scale_factor as choose_scale_factor,
)
from nwkit.calibration_scale import (
    cmd_factor as cmd_factor,
)
from nwkit.calibration_scale import (
    cmd_scale_calibrations as cmd_scale_calibrations,
)
from nwkit.calibration_scale import (
    cmd_scale_ctl_rootage as cmd_scale_ctl_rootage,
)
from nwkit.calibration_scale import (
    format_decimal as format_decimal,
)
from nwkit.calibration_scale import (
    main as main,
)
from nwkit.calibration_scale import (
    parse_decimal as parse_decimal,
)
from nwkit.calibration_scale import (
    parse_scale as parse_scale,
)
from nwkit.calibration_scale import (
    read_text as read_text,
)
from nwkit.calibration_scale import (
    scale_calibration_labels as scale_calibration_labels,
)
from nwkit.calibration_scale import (
    scale_ctl_rootage_text as scale_ctl_rootage_text,
)
from nwkit.calibration_scale import (
    scale_decimal as scale_decimal,
)
from nwkit.calibration_scale import (
    scale_number_text as scale_number_text,
)
from nwkit.calibration_scale import (
    scale_rootage_line as scale_rootage_line,
)
from nwkit.calibration_scale import (
    split_field_spacing as split_field_spacing,
)
from nwkit.calibration_scale import (
    write_text as write_text,
)

if __name__ == "__main__":
    raise SystemExit(main())

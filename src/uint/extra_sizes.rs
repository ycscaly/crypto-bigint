//! Support for additional integer sizes beyond the core set which is defined
//! in the toplevel module.
//!
//! These are feature-gated to keep compile times down for applications which
//! do not need them.
// TODO(tarcieri): switch to a fully const generic implementation using `generic_const_exprs`

use super::*;

impl_uint_aliases! {
    (U704, 704, "704-bit"),
    (U832, 832, "832-bit"),
    (U960, 960, "960-bit"),
    (U1088, 1088, "1088-bit"),
    (U1152, 1152, "1152-bit"),
    (U1216, 1216, "1216-bit"),
    (U1344, 1344, "1344-bit"),
    (U1408, 1408, "1408-bit"),
    (U1472, 1472, "1472-bit"),
    (U1600, 1600, "1600-bit"),
    (U1664, 1664, "1664-bit"),
    (U1728, 1728, "1728-bit"),
    (U1856, 1856, "1856-bit"),
    (U1920, 1920, "1920-bit"),
    (U1984, 1984, "1984-bit"),
    (U2112, 2112, "2112-bit"),
    (U2176, 2176, "2176-bit"),
    (U2240, 2240, "2240-bit"),
    (U2304, 2304, "2304-bit"),
    (U2368, 2368, "2368-bit"),
    (U2432, 2432, "2432-bit"),
    (U2496, 2496, "2496-bit"),
    (U2560, 2560, "2560-bit"),
    (U2624, 2624, "2624-bit"),
    (U2688, 2688, "2688-bit"),
    (U2752, 2752, "2752-bit"),
    (U2816, 2816, "2816-bit"),
    (U2880, 2880, "2880-bit"),
    (U2944, 2944, "2944-bit"),
    (U3008, 3008, "3008-bit"),
    (U3136, 3136, "3136-bit"),
    (U3200, 3200, "3200-bit"),
    (U3264, 3264, "3264-bit"),
    (U3328, 3328, "3328-bit"),
    (U3392, 3392, "3392-bit"),
    (U3456, 3456, "3456-bit"),
    (U3520, 3520, "3520-bit"),
    (U3648, 3648, "3648-bit"),
    (U3712, 3712, "3712-bit"),
    (U3776, 3776, "3776-bit"),
    (U3840, 3840, "3840-bit"),
    (U3904, 3904, "3904-bit"),
    (U3968, 3968, "3968-bit"),
    (U4032, 4032, "4032-bit"),
    (U4160, 4160, "4160-bit"),
    (U4288, 4288, "4288-bit"),
    (U4416, 4416, "4416-bit"),
    (U4480, 4480, "4480-bit"),
    (U4544, 4544, "4544-bit"),
    (U4608, 4608, "4608-bit"),
    (U4672, 4672, "4672-bit"),
    (U4736, 4736, "4736-bit"),
    (U4800, 4800, "4800-bit"),
    (U4864, 4864, "4864-bit"),
    (U4928, 4928, "4928-bit"),
    (U4992, 4992, "4992-bit"),
    (U5056, 5056, "5056-bit"),
    (U5120, 5120, "5120-bit"),
    (U5184, 5184, "5184-bit"),
    (U5248, 5248, "5248-bit"),
    (U5312, 5312, "5312-bit"),
    (U5376, 5376, "5376-bit"),
    (U5440, 5440, "5440-bit"),
    (U5504, 5504, "5504-bit"),
    (U5568, 5568, "5568-bit"),
    (U5632, 5632, "5632-bit"),
    (U5696, 5696, "5696-bit"),
    (U5760, 5760, "5760-bit"),
    (U5824, 5824, "5824-bit"),
    (U5888, 5888, "5888-bit"),
    (U5952, 5952, "5952-bit"),
    (U6016, 6016, "6016-bit"),
    (U6080, 6080, "6080-bit"),
    (U6208, 6208, "6208-bit"),
    (U6272, 6272, "6272-bit"),
    (U6336, 6336, "6336-bit"),
    (U6400, 6400, "6400-bit"),
    (U6464, 6464, "6464-bit"),
    (U6528, 6528, "6528-bit"),
    (U6592, 6592, "6592-bit"),
    (U6656, 6656, "6656-bit"),
    (U6720, 6720, "6720-bit"),
    (U6784, 6784, "6784-bit"),
    (U6848, 6848, "6848-bit"),
    (U6912, 6912, "6912-bit"),
    (U6976, 6976, "6976-bit"),
    (U7040, 7040, "7040-bit"),
    (U7104, 7104, "7104-bit"),
    (U7168, 7168, "7168-bit"),
    (U7232, 7232, "7232-bit"),
    (U7296, 7296, "7296-bit"),
    (U7360, 7360, "7360-bit"),
    (U7424, 7424, "7424-bit"),
    (U7488, 7488, "7488-bit"),
    (U7552, 7552, "7552-bit"),
    (U7616, 7616, "7616-bit"),
    (U7680, 7680, "7680-bit"),
    (U7744, 7744, "7744-bit"),
    (U7808, 7808, "7808-bit"),
    (U7872, 7872, "7872-bit"),
    (U7936, 7936, "7936-bit"),
    (U8000, 8000, "8000-bit"),
    (U8064, 8064, "8064-bit"),
    (U8128, 8128, "8128-bit")
}

impl_concat! {
    (U576, 576),
    (U704, 704),
    (U832, 832),
    (U960, 960),
    (U1088, 1088),
    (U1152, 1152),
    (U1216, 1216),
    (U1280, 1280),
    (U1344, 1344),
    (U1408, 1408),
    (U1472, 1472),
    (U1600, 1600),
    (U1664, 1664),
    (U1728, 1728),
    (U1856, 1856),
    (U1920, 1920),
    (U1984, 1984),
    (U2112, 2112),
    (U2176, 2176),
    (U2240, 2240),
    (U2304, 2304),
    (U2368, 2368),
    (U2432, 2432),
    (U2496, 2496),
    (U2560, 2560),
    (U2624, 2624),
    (U2688, 2688),
    (U2752, 2752),
    (U2816, 2816),
    (U2880, 2880),
    (U2944, 2944),
    (U3008, 3008),
    (U3136, 3136),
    (U3200, 3200),
    (U3264, 3264),
    (U3328, 3328),
    (U3392, 3392),
    (U3456, 3456),
    (U3520, 3520),
    (U3584, 3584),
    (U3648, 3648),
    (U3712, 3712),
    (U3776, 3776),
    (U3840, 3840),
    (U3904, 3904),
    (U3968, 3968),
    (U4032, 4032),
    (U4160, 4160),
    (U4288, 4288),
    (U4416, 4416),
    (U4480, 4480),
    (U4544, 4544),
    (U4608, 4608),
    (U4672, 4672),
    (U4736, 4736),
    (U4800, 4800),
    (U4864, 4864),
    (U4928, 4928),
    (U4992, 4992),
    (U5056, 5056),
    (U5120, 5120),
    (U5184, 5184),
    (U5248, 5248),
    (U5312, 5312),
    (U5376, 5376),
    (U5440, 5440),
    (U5504, 5504),
    (U5568, 5568),
    (U5632, 5632),
    (U5696, 5696),
    (U5760, 5760),
    (U5824, 5824),
    (U5888, 5888),
    (U5952, 5952),
    (U6016, 6016),
    (U6080, 6080),
    (U6144, 6144),
    (U6208, 6208),
    (U6272, 6272),
    (U6336, 6336),
    (U6400, 6400),
    (U6464, 6464),
    (U6528, 6528),
    (U6592, 6592),
    (U6656, 6656),
    (U6720, 6720),
    (U6784, 6784),
    (U6848, 6848),
    (U6912, 6912),
    (U6976, 6976),
    (U7040, 7040),
    (U7104, 7104),
    (U7168, 7168),
    (U7232, 7232),
    (U7296, 7296),
    (U7360, 7360),
    (U7424, 7424),
    (U7488, 7488),
    (U7552, 7552),
    (U7616, 7616),
    (U7680, 7680),
    (U7744, 7744),
    (U7808, 7808),
    (U7872, 7872),
    (U7936, 7936),
    (U8000, 8000),
    (U8064, 8064),
    (U8128, 8128)
}

impl_split! {
    (U1152, 1152),
    (U1408, 1408),
    (U1664, 1664),
    (U1920, 1920),
    (U2176, 2176),
    (U2304, 2304),
    (U2432, 2432),
    (U2560, 2560),
    (U2688, 2688),
    (U2816, 2816),
    (U2944, 2944),
    (U3200, 3200),
    (U3328, 3328),
    (U3456, 3456),
    (U3712, 3712),
    (U3840, 3840),
    (U3968, 3968),
    (U4480, 4480),
    (U4608, 4608),
    (U4736, 4736),
    (U4864, 4864),
    (U4992, 4992),
    (U5120, 5120),
    (U5248, 5248),
    (U5376, 5376),
    (U5504, 5504),
    (U5632, 5632),
    (U5760, 5760),
    (U5888, 5888),
    (U6016, 6016),
    (U6272, 6272),
    (U6400, 6400),
    (U6528, 6528),
    (U6656, 6656),
    (U6784, 6784),
    (U6912, 6912),
    (U7040, 7040),
    (U7168, 7168),
    (U7296, 7296),
    (U7424, 7424),
    (U7552, 7552),
    (U7680, 7680),
    (U7808, 7808),
    (U7936, 7936),
    (U8064, 8064)
}

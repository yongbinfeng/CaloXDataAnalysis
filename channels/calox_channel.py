from channels.gain_validator import enforce_gain
import numpy as np
import copy
import logging

# Set up logging for a cleaner output
logging.basicConfig(level=logging.WARNING,
                    format='\033[91m%(levelname)s: %(message)s\033[0m')


class CaloXChannel(object):
    """
    A class to represent a channel in the FERS system.
    """

    def __init__(self, i_tower_x: float, i_tower_y: float, isCer: bool, isQuartz: bool = False, is6mm: bool = True):
        self.i_tower_x = i_tower_x
        self.i_tower_y = i_tower_y
        self.isCer = isCer
        self.isQuartz = isQuartz
        self.is6mm = is6mm

    def __str__(self):
        # Use dynamic class name + vars
        attrs = ', '.join(f"{k}={v}" for k, v in vars(self).items())
        return f"{self.__class__.__name__}({attrs})"

    def __copy__(self):
        # Works for subclasses automatically
        cls = self.__class__
        new_obj = cls.__new__(cls)
        new_obj.__dict__.update(self.__dict__)
        return new_obj

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.__dict__ == other.__dict__

    def isCer(self):
        return self.isCer

    def isSci(self):
        return not self.isCer

    def isQuartz(self):
        return self.isQuartz and self.isCer


class FERSChannel(CaloXChannel):
    def __init__(self, i_tower_x: float, i_tower_y: float, isCer: bool, channel_no: int, board_no: int, isQuartz: bool = False, is6mm: bool = True):
        super().__init__(i_tower_x, i_tower_y, isCer, isQuartz, is6mm)
        self.channel_no = channel_no  # FERS channel number
        self.board_no = board_no  # FERS board number

    @enforce_gain
    def get_channel_name(self, gain="HG", pdsub=False, calib=False):
        channelName = f"FERS_Board{self.board_no}_energy{gain}_{self.channel_no}"
        if pdsub and gain != "Mix":
            # mix channels always have pedestal subtracted
            channelName += "_pdsub"
        if calib:
            channelName += "_calib"
        return channelName

    def get_tower_pos_name(self, isX=True):
        return f"FERS_Board{self.board_no}_{self.channel_no}_iTower{'X' if isX else 'Y'}"

    def get_real_pos_name(self, isX=True):
        return f"FERS_Board{self.board_no}_{self.channel_no}_Real{'X' if isX else 'Y'}"


class DRSChannel(CaloXChannel):
    def __init__(self, i_tower_x: float, i_tower_y: float, isCer: bool,
                 channel_no: int, group_no: int, board_no: int,
                 is_amplified: bool = False, is6mm: bool = True, isQuartz: bool = False):
        super().__init__(i_tower_x, i_tower_y, isCer, isQuartz, is6mm)
        self.channel_no = channel_no
        self.group_no = group_no
        self.board_no = board_no
        self.is_amplified = is_amplified

    def get_channel_name(self, blsub=False):
        channelName = f"DRS_Board{self.board_no}_Group{self.group_no}_Channel{self.channel_no}"
        if blsub:
            channelName += "_blsub"
        return channelName

    def get_channel_sum_name(self):
        return self.get_channel_name(blsub=True) + "_sum"

    def get_channel_peak_name(self):
        return self.get_channel_name(blsub=True) + "_peak"

    def get_channel_peak_ts_name(self):
        return self.get_channel_name(blsub=True) + "_peakTS"


class Board(object):
    """
    A class to represent a generic board.
    """

    def __init__(self, board_no):
        self.board_no = board_no  # Board number
        self.channels = []  # List of channels on the board

    def __str__(self):
        s = f"Board Number: {self.board_no}\n"
        for channel in self.channels:
            s += str(channel) + "\n"
        return s.strip()

    def __iter__(self):
        for channel in self.channels:
            yield channel

    def __copy__(self):
        cls = self.__class__
        new_board = cls.__new__(cls)
        new_board.__dict__.update(self.__dict__)
        new_board.channels = [copy.copy(channel) for channel in self.channels]
        return new_board

    def copy(self, board_no):
        # Create a new instance of the same class (e.g., FERSBoard or DRSBoard)
        new_board = self.__class__.__new__(self.__class__)
        # Copy all attributes from the old object to the new object
        new_board.__dict__.update(self.__dict__)
        new_board.board_no = board_no

        new_board.channels = [copy.deepcopy(
            channel) for channel in self.channels]
        for channel in new_board.channels:
            if hasattr(channel, 'board_no'):
                channel.board_no = board_no
        return new_board

    def get_channel_by_tower(self, i_tower_x, i_tower_y, isCer=False):
        """
        Get a channel by its tower coordinates (iTowerX, iTowerY).
        Returns the first matching channel or None if not found.
        """
        for channel in self.channels:
            if channel.i_tower_x == i_tower_x and channel.i_tower_y == i_tower_y and channel.isCer == isCer:
                return channel
        # logging.warning(
        #    f"Channel not found for tower ({iTowerX}, {iTowerY}) on board {self.boardNo} for {'CER' if isCer else 'SCI'} channel.")
        return None

    def get_list_of_channels(self, isCer=None):
        """
        isCer can be True, False, or None (for all channels).
        """
        if isCer is None:
            return self.channels
        return [channel for channel in self.channels if channel.isCer == isCer]

    def get_list_of_towers(self):
        """
        Get a list of unique towers on the board.
        Returns a list of tuples (iTowerX, iTowerY).
        """
        return list(set((channel.i_tower_x, channel.i_tower_y) for channel in self.channels))

    def map_cer_to_sci(self):
        """
        Map CER channels to SCI channels.
        Returns a dictionary mapping CER channels to SCI channels.
        """
        cer_channels = self.get_list_of_channels(isCer=True)
        sci_channels = self.get_list_of_channels(isCer=False)
        mapping = {}
        for cer_channel in cer_channels:
            for sci_channel in sci_channels:
                if cer_channel.i_tower_x == sci_channel.i_tower_x and \
                   cer_channel.i_tower_y == sci_channel.i_tower_y:
                    mapping[cer_channel] = sci_channel
                    break
        assert len(mapping) == len(cer_channels), \
            "Mapping failed: not all CER channels have a corresponding SCI channel."
        assert len(set(mapping.values())) == len(mapping), \
            "Mapping failed: some SCI channels are mapped to multiple CER channels."
        return mapping

    def shift_channels(self, i_shift_x, i_shift_y):
        """
        Shift the channels on the board by (iShiftX, iShiftY).
        """
        for channel in self.channels:
            channel.i_tower_x += i_shift_x
            channel.i_tower_y += i_shift_y

    def move_to(self, i_tower_x, i_tower_y):
        """
        Move the channels on the board to a new position.
        Use the first channel as reference.
        """
        if not self.channels:
            logging.warning(f"No channels to move on board {self.board_no}.")
            return
        shift_x = i_tower_x - self.channels[0].i_tower_x
        shift_y = i_tower_y - self.channels[0].i_tower_y
        self.shift_channels(shift_x, shift_y)


class FERSBoard(Board):
    """
    A class to represent a FERS board.
    """

    def __init__(self, board_no, is6mm=False):
        """
        Initialize a FERS board with a specific board number and channel configuration.
        :param boardNo: The board number.
        :param is6mm: Boolean indicating if the board is 6mm (True) or 3mm (False).
        """
        super().__init__(board_no)
        self.is6mm = is6mm  # Flag for 6mm or 3mm board
        # channels is a 2D list of FERSChannel objects
        self.channels = build_fers_base(is6mm=is6mm, board_no=board_no)

    def is_6mm_size(self):
        return self.is6mm

    def is_3mm_size(self):
        return not self.is_6mm_size()

    @enforce_gain
    def get_energy_max_name(self, gain="HG", isCer=True):
        cat = "Cer" if isCer else "Sci"
        # name is based on the channel name
        # "FERS_Board{self.boardNo}_energyHG_{self.channelNo}"
        return f"FERS_Board{self.board_no}_energy{gain}_{cat}_max"

    @enforce_gain
    def get_energy_sum_name(self, gain="HG", isCer=True, pdsub=False, calib=False):
        cat = "Cer" if isCer else "Sci"
        sumname = f"FERS_Board{self.board_no}_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            sumname += "_pdsub"
        if calib:
            sumname += "_calib"
        sumname += "_sum"
        return sumname

    @enforce_gain
    def get_energy_weighted_center_name(self, gain="HG", isCer=True, pdsub=False, calib=False, isX=True):
        cat = "Cer" if isCer else "Sci"
        centername = f"FERS_Board{self.board_no}_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            centername += "_pdsub"
        if calib:
            centername += "_calib"
        centername += "_weighted_center"
        centername += "_X" if isX else "_Y"
        return centername

    def get_sipm_hv_name(self):
        return f"FERS_Board{self.board_no}_SipmHV"

    def get_sipm_i_name(self):
        return f"FERS_Board{self.board_no}_SipmI"

    def get_temp_det_name(self):
        return f"FERS_Board{self.board_no}_TempDET"

    def get_temp_fpga_name(self):
        return f"FERS_Board{self.board_no}_TempFPGA"


class DRSBoard(Board):
    """
    A class to represent a DRS board.
    channels: flat list[DRSChannel]
    """

    def __init__(self, board_no, is6mm=True, channels=None):
        """
        Initialize a DRS board with a specific board number and channel configuration.
        :param boardNo: The board number.
        """
        super().__init__(board_no)
        # channels is a 2D list of DRSChannel objects
        if channels is not None:
            self.channels = channels
        else:
            self.channels = build_drs_base(is6mm=is6mm, board_no=board_no)

    def get_channel_by_group_channel(self, group_no, chan_no):
        """
        Get a channel by group number and channel number.
        Returns the first matching channel or None if not found.
        """
        for channel in self.channels:
            if channel.group_no == group_no and channel.channel_no == chan_no:
                return channel
        logging.warning(
            f"Channel Group{group_no} Channel{chan_no} not found on board {self.board_no}.")
        return None

    def remove_channel_by_group_channel(self, group_no, chan_no):
        """
        Remove a channel from the board by group number and channel number.
        """
        channel_to_remove = self.get_channel_by_group_channel(group_no, chan_no)
        if channel_to_remove:
            self.channels.remove(channel_to_remove)
        else:
            logging.warning(
                f"Channel Group{group_no} Channel{chan_no} not found on board {self.board_no}.")


class FERSBoards(dict):
    """
    A dictionary to hold multiple FERS boards.
    Key: boardNo, Value: FERSBoard object
    """

    def __setitem__(self, key, value):
        if not isinstance(value, FERSBoard):
            raise ValueError("Value must be a FERSBoard instance.")
        super().__setitem__(key, value)

    @enforce_gain
    def get_energy_max_name(self, gain="HG", isCer=True):
        cat = "Cer" if isCer else "Sci"
        return f"FERS_energy{gain}_{cat}_max"

    @enforce_gain
    def get_energy_sum_name(self, gain="HG", isCer=True, pdsub=False, calib=False):
        cat = "Cer" if isCer else "Sci"
        sumname = f"FERS_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            sumname += "_pdsub"
        if calib:
            sumname += "_calib"
        sumname += "_sum"
        return sumname

    @enforce_gain
    def get_energy_weighted_center_name(self, gain="HG", isCer=True, pdsub=False, calib=False, isX=True):
        cat = "Cer" if isCer else "Sci"
        centername = f"FERS_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            centername += "_pdsub"
        if calib:
            centername += "_calib"
        centername += "_weighted_center"
        centername += "_X" if isX else "_Y"
        return centername


# physical channels to the readout channel names
A5202_map = np.swapaxes(
    np.array([
        [0, 2, 1, 3],
        [4, 6, 5, 7],
        [10, 8, 11, 9],
        [14, 12, 15, 13],
        [16, 18, 17, 19],
        [20, 22, 21, 23],
        [26, 24, 27, 25],
        [30, 28, 31, 29],
        [32, 34, 33, 35],
        [36, 38, 37, 39],
        [42, 40, 43, 41],
        [46, 44, 47, 45],
        [48, 50, 49, 51],
        [52, 54, 53, 55],
        [58, 56, 59, 57],
        [62, 60, 63, 61]
    ], dtype=int), 0, 1)

A5205_map_3mm = np.swapaxes(
    np.array([
        [2, 0, 1, 3],
        [6, 4, 5, 7],
        [10, 8, 9, 11],
        [14, 12, 13, 15],
        [16, 18, 19, 17],
        [20, 22, 23, 21],
        [24, 26, 27, 25],
        [28, 30, 31, 29],
        [32, 34, 35, 33],
        [36, 38, 39, 37],
        [40, 42, 43, 41],
        [44, 46, 47, 45],
        [48, 50, 51, 49],
        [52, 54, 55, 53],
        [56, 58, 59, 57],
        [60, 62, 63, 61]
    ], dtype=int), 0, 1)


drs_map = np.swapaxes(
    np.array([
        [3, 2, 1, 0],
        [7, 6, 5, 4],
        [11, 10, 9, 8],
        [15, 14, 13, 12],
        [19, 18, 17, 16],
        [23, 22, 21, 20],
        [27, 26, 25, 24],
        [31, 30, 29, 28],
    ], dtype=int), 0, 1)


def build_fers_base(is6mm=False, board_no=0):
    channels_FERS = []
    for ix in range(0, 4):
        for iy in range(0, 16):
            if is6mm:
                channel_no = A5202_map[ix, iy]
                isCer = (iy % 2 == 0)
                channel = FERSChannel(ix, -int(iy/2), isCer,
                                      channel_no, board_no)
                channels_FERS.append(channel)
            else:
                channel_no = A5205_map_3mm[ix, iy]
                isCer = (ix % 2 == 0)
                # 3mm has higher granularity in global Y
                channel = FERSChannel(
                    int(ix/2), -float(iy)/4, isCer, channel_no, board_no)
                channels_FERS.append(channel)
    return channels_FERS


def build_drs_base(is6mm=True, board_no=0):
    # right now only have DRS on 6mm
    channels_DRS = []
    for ix in range(0, 4):
        for iy in range(0, 8):
            channel_no = drs_map[ix, iy]
            # DRS channels are grouped in groups of 8
            # groupNo is the group number (0-3)
            # chanNo is the channel number within the group (0-7)
            group_no = (channel_no // 8)
            chan_no = (channel_no % 8)
            if is6mm:
                isCer = (iy % 2 == 0)
                channel = DRSChannel(ix, -int(iy/2), isCer,
                                     chan_no, group_no, board_no, is6mm=is6mm)
            else:
                # this is INCONSISTENT with the FERS 3mm channels
                # need to CHECK
                isCer = (ix % 2 == 1)
                # 3mm has higher granularity in global Y
                channel = DRSChannel(
                    int(ix/2), -float(iy)/4, isCer, chan_no, group_no, board_no, is6mm=False, is_amplified=True)
            channels_DRS.append(channel)
    return channels_DRS


if __name__ == "__main__":
    # Example usage
    print("Building FERS board 6mm:")
    fers_board = FERSBoard(board_no=1, is6mm=True)
    print(fers_board)

    print("\nBuilding FERS board 3mm:")
    fers_board_3mm = FERSBoard(board_no=2, is6mm=False)
    print(fers_board_3mm)

    print("\nBuilding DRS board:")
    drs_board = DRSBoard(board_no=1)
    print(drs_board)

    print("\nDRS Board List of Towers: ", drs_board.get_list_of_towers())

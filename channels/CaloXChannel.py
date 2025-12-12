from channels.gainvalidator import enforce_gain
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

    def __init__(self, iTowerX: float, iTowerY: float, isCer: bool, isQuartz: bool = False, is6mm: bool = True):
        self.iTowerX = iTowerX
        self.iTowerY = iTowerY
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
    def __init__(self, iTowerX: float, iTowerY: float, isCer: bool, channelNo: int, boardNo: int, isQuartz: bool = False, is6mm: bool = True):
        super().__init__(iTowerX, iTowerY, isCer, isQuartz, is6mm)
        self.channelNo = channelNo  # FERS channel number
        self.boardNo = boardNo  # FERS board number

    @enforce_gain
    def GetChannelName(self, gain="HG", pdsub=False, calib=False):
        channelName = f"FERS_Board{self.boardNo}_energy{gain}_{self.channelNo}"
        if pdsub and gain != "Mix":
            # mix channels always have pedestal subtracted
            channelName += "_pdsub"
        if calib:
            channelName += "_calib"
        return channelName

    def GetTowerPosName(self, isX=True):
        return f"FERS_Board{self.boardNo}_{self.channelNo}_iTower{'X' if isX else 'Y'}"

    def GetRealPosName(self, isX=True):
        return f"FERS_Board{self.boardNo}_{self.channelNo}_Real{'X' if isX else 'Y'}"


class DRSChannel(CaloXChannel):
    def __init__(self, iTowerX: float, iTowerY: float, isCer: bool,
                 channelNo: int, groupNo: int, boardNo: int,
                 isAmplified: bool = False, is6mm: bool = True, isQuartz: bool = False):
        super().__init__(iTowerX, iTowerY, isCer, isQuartz, is6mm)
        self.channelNo = channelNo
        self.groupNo = groupNo
        self.boardNo = boardNo
        self.isAmplified = isAmplified

    def GetChannelName(self, blsub=False):
        channelName = f"DRS_Board{self.boardNo}_Group{self.groupNo}_Channel{self.channelNo}"
        if blsub:
            channelName += "_blsub"
        return channelName

    def GetChannelSumName(self):
        return self.GetChannelName(blsub=True) + "_sum"

    def GetChannelPeakName(self):
        return self.GetChannelName(blsub=True) + "_peak"

    def GetChannelPeakTSName(self):
        return self.GetChannelName(blsub=True) + "_peakTS"


class Board(object):
    """
    A class to represent a generic board.
    """

    def __init__(self, boardNo):
        self.boardNo = boardNo  # Board number
        self.channels = []  # List of channels on the board

    def __str__(self):
        s = f"Board Number: {self.boardNo}\n"
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

    def copy(self, boardNo):
        # Create a new instance of the same class (e.g., FERSBoard or DRSBoard)
        new_board = self.__class__.__new__(self.__class__)
        # Copy all attributes from the old object to the new object
        new_board.__dict__.update(self.__dict__)
        new_board.boardNo = boardNo

        new_board.channels = [copy.deepcopy(
            channel) for channel in self.channels]
        for channel in new_board.channels:
            if hasattr(channel, 'boardNo'):
                channel.boardNo = boardNo
        return new_board

    def GetChannelByTower(self, iTowerX, iTowerY, isCer=False):
        """
        Get a channel by its tower coordinates (iTowerX, iTowerY).
        Returns the first matching channel or None if not found.
        """
        for channel in self.channels:
            if channel.iTowerX == iTowerX and channel.iTowerY == iTowerY and channel.isCer == isCer:
                return channel
        # logging.warning(
        #    f"Channel not found for tower ({iTowerX}, {iTowerY}) on board {self.boardNo} for {'CER' if isCer else 'SCI'} channel.")
        return None

    def GetListOfChannels(self, isCer=None):
        """
        isCer can be True, False, or None (for all channels).
        """
        if isCer is None:
            return self.channels
        return [channel for channel in self.channels if channel.isCer == isCer]

    def GetListOfTowers(self):
        """
        Get a list of unique towers on the board.
        Returns a list of tuples (iTowerX, iTowerY).
        """
        return list(set((channel.iTowerX, channel.iTowerY) for channel in self.channels))

    def MapCerToSci(self):
        """
        Map CER channels to SCI channels.
        Returns a dictionary mapping CER channels to SCI channels.
        """
        cer_channels = self.GetListOfChannels(isCer=True)
        sci_channels = self.GetListOfChannels(isCer=False)
        mapping = {}
        for cer_channel in cer_channels:
            for sci_channel in sci_channels:
                if cer_channel.iTowerX == sci_channel.iTowerX and \
                   cer_channel.iTowerY == sci_channel.iTowerY:
                    mapping[cer_channel] = sci_channel
                    break
        assert len(mapping) == len(cer_channels), \
            "Mapping failed: not all CER channels have a corresponding SCI channel."
        assert len(set(mapping.values())) == len(mapping), \
            "Mapping failed: some SCI channels are mapped to multiple CER channels."
        return mapping

    def ShiftChannels(self, iShiftX, iShiftY):
        """
        Shift the channels on the board by (iShiftX, iShiftY).
        """
        for channel in self.channels:
            channel.iTowerX += iShiftX
            channel.iTowerY += iShiftY

    def MoveTo(self, iTowerX, iTowerY):
        """
        Move the channels on the board to a new position.
        Use the first channel as reference.
        """
        if not self.channels:
            logging.warning(f"No channels to move on board {self.boardNo}.")
            return
        shiftX = iTowerX - self.channels[0].iTowerX
        shiftY = iTowerY - self.channels[0].iTowerY
        self.ShiftChannels(shiftX, shiftY)


class FERSBoard(Board):
    """
    A class to represent a FERS board.
    """

    def __init__(self, boardNo, is6mm=False):
        """
        Initialize a FERS board with a specific board number and channel configuration.
        :param boardNo: The board number.
        :param is6mm: Boolean indicating if the board is 6mm (True) or 3mm (False).
        """
        super().__init__(boardNo)
        self.is6mm = is6mm  # Flag for 6mm or 3mm board
        # channels is a 2D list of FERSChannel objects
        self.channels = buildFERSBase(is6mm=is6mm, boardNo=boardNo)

    def Is6mm(self):
        return self.is6mm

    def Is3mm(self):
        return not self.Is6mm()

    @enforce_gain
    def GetEnergyMaxName(self, gain="HG", isCer=True):
        cat = "Cer" if isCer else "Sci"
        # name is based on the channel name
        # "FERS_Board{self.boardNo}_energyHG_{self.channelNo}"
        return f"FERS_Board{self.boardNo}_energy{gain}_{cat}_max"

    @enforce_gain
    def GetEnergySumName(self, gain="HG", isCer=True, pdsub=False, calib=False):
        cat = "Cer" if isCer else "Sci"
        sumname = f"FERS_Board{self.boardNo}_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            sumname += "_pdsub"
        if calib:
            sumname += "_calib"
        sumname += "_sum"
        return sumname

    @enforce_gain
    def GetEnergyWeightedCenterName(self, gain="HG", isCer=True, pdsub=False, calib=False, isX=True):
        cat = "Cer" if isCer else "Sci"
        centername = f"FERS_Board{self.boardNo}_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            centername += "_pdsub"
        if calib:
            centername += "_calib"
        centername += "_weighted_center"
        centername += "_X" if isX else "_Y"
        return centername

    def GetSipmHVName(self):
        return f"FERS_Board{self.boardNo}_SipmHV"

    def GetSipmIName(self):
        return f"FERS_Board{self.boardNo}_SipmI"

    def GetTempDETName(self):
        return f"FERS_Board{self.boardNo}_TempDET"

    def GetTempFPGAName(self):
        return f"FERS_Board{self.boardNo}_TempFPGA"


class DRSBoard(Board):
    """
    A class to represent a DRS board.
    channels: flat list[DRSChannel]
    """

    def __init__(self, boardNo, is6mm=True, channels=None):
        """
        Initialize a DRS board with a specific board number and channel configuration.
        :param boardNo: The board number.
        """
        super().__init__(boardNo)
        # channels is a 2D list of DRSChannel objects
        if channels is not None:
            self.channels = channels
        else:
            self.channels = buildDRSBase(is6mm=is6mm, boardNo=boardNo)

    def GetChannelByGroupChannel(self, groupNo, chanNo):
        """
        Get a channel by group number and channel number.
        Returns the first matching channel or None if not found.
        """
        for channel in self.channels:
            if channel.groupNo == groupNo and channel.channelNo == chanNo:
                return channel
        logging.warning(
            f"Channel Group{groupNo} Channel{chanNo} not found on board {self.boardNo}.")
        return None

    def RemoveChannelByGroupChannel(self, groupNo, chanNo):
        """
        Remove a channel from the board by group number and channel number.
        """
        channel_to_remove = self.GetChannelByGroupChannel(groupNo, chanNo)
        if channel_to_remove:
            self.channels.remove(channel_to_remove)
        else:
            logging.warning(
                f"Channel Group{groupNo} Channel{chanNo} not found on board {self.boardNo}.")


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
    def GetEnergyMaxName(self, gain="HG", isCer=True):
        cat = "Cer" if isCer else "Sci"
        return f"FERS_energy{gain}_{cat}_max"

    @enforce_gain
    def GetEnergySumName(self, gain="HG", isCer=True, pdsub=False, calib=False):
        cat = "Cer" if isCer else "Sci"
        sumname = f"FERS_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            sumname += "_pdsub"
        if calib:
            sumname += "_calib"
        sumname += "_sum"
        return sumname

    @enforce_gain
    def GetEnergySumRatioName(self, gain="HG", pdsub=False, calib=False):
        cat = "CerOverSci"
        ratioName = f"FERS_energy{gain}_{cat}"
        if pdsub and gain != "Mix":
            ratioName += "_pdsub"
        if calib:
            ratioName += "_calib"
        ratioName += "_sum"
        return ratioName

    @enforce_gain
    def GetEnergyWeightedCenterName(self, gain="HG", isCer=True, pdsub=False, calib=False, isX=True):
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


def buildFERSBase(is6mm=False, boardNo=0):
    channels_FERS = []
    for ix in range(0, 4):
        for iy in range(0, 16):
            if is6mm:
                channelNo = A5202_map[ix, iy]
                isCer = (iy % 2 == 0)
                channel = FERSChannel(ix, -int(iy/2), isCer,
                                      channelNo, boardNo)
                channels_FERS.append(channel)
            else:
                channelNo = A5205_map_3mm[ix, iy]
                isCer = (ix % 2 == 0)
                # 3mm has higher granularity in global Y
                channel = FERSChannel(
                    int(ix/2), -float(iy)/4, isCer, channelNo, boardNo)
                channels_FERS.append(channel)
    return channels_FERS


def buildDRSBase(is6mm=True, boardNo=0):
    # right now only have DRS on 6mm
    channels_DRS = []
    for ix in range(0, 4):
        for iy in range(0, 8):
            channelNo = drs_map[ix, iy]
            # DRS channels are grouped in groups of 8
            # groupNo is the group number (0-3)
            # chanNo is the channel number within the group (0-7)
            groupNo = (channelNo // 8)
            chanNo = (channelNo % 8)
            if is6mm:
                isCer = (iy % 2 == 0)
                channel = DRSChannel(ix, -int(iy/2), isCer,
                                     chanNo, groupNo, boardNo, is6mm=is6mm)
            else:
                # this is INCONSISTENT with the FERS 3mm channels
                # need to CHECK
                isCer = (ix % 2 == 1)
                # 3mm has higher granularity in global Y
                channel = DRSChannel(
                    int(ix/2), -float(iy)/4, isCer, chanNo, groupNo, boardNo, is6mm=False, isAmplified=True)
            channels_DRS.append(channel)
    return channels_DRS


if __name__ == "__main__":
    # Example usage
    print("Building FERS board 6mm:")
    fers_board = FERSBoard(boardNo=1, is6mm=True)
    print(fers_board)

    print("\nBuilding FERS board 3mm:")
    fers_board_3mm = FERSBoard(boardNo=2, is6mm=False)
    print(fers_board_3mm)

    print("\nBuilding DRS board:")
    drs_board = DRSBoard(boardNo=1)
    print(drs_board)

    print("\nDRS Board List of Towers: ", drs_board.GetListOfTowers())

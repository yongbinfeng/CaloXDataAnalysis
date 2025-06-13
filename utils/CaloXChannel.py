import numpy as np


class CaloXChannel(object):
    """
    A class to represent a channel in the FERS system.
    """

    def __init__(self, iTowerX: int, iTowerY: int, iBoardX: int, iBoardY: int, isCer: bool):
        self.iTowerX = iTowerX
        self.iTowerY = iTowerY
        self.iBoardX = iBoardX
        self.iBoardY = iBoardY
        self.isCer = isCer

    def __str__(self):
        return (f"CaloXChannel(iTowerX={self.iTowerX}, "
                f"iTowerY={self.iTowerY}, iBoardX={self.iBoardX}, "
                f"iBoardY={self.iBoardY}, isCer={self.isCer}")

    def __copy__(self):
        """
        Create a shallow copy of the channel.
        """
        return CaloXChannel(
            self.iTowerX, self.iTowerY, self.iBoardX, self.iBoardY,
            self.isCer
        )

    def isCer(self):
        return self.isCer

    def isSci(self):
        return not self.isCer


class FERSChannel(CaloXChannel):
    def __init__(self, iTowerX: int, iTowerY: int, iBoardX: int, iBoardY: int, isCer: bool, channelNo: int, boardNo: int):
        super().__init__(iTowerX, iTowerY, iBoardX, iBoardY, isCer)
        self.channelNo = channelNo  # FERS channel number
        self.boardNo = boardNo  # FERS board number

    def __str__(self):
        base_str = super().__str__()
        base_str += f", channelNo={self.channelNo}, boardNo={self.boardNo})"
        base_str = base_str.replace("CaloXChannel", "FERSChannel")
        return base_str

    def __copy__(self):
        """
        Create a shallow copy of the FERS channel.
        """
        return FERSChannel(
            self.iTowerX, self.iTowerY, self.iBoardX, self.iBoardY,
            self.isCer, self.channelNo, self.boardNo
        )

    def GetHGChannelName(self):
        return f"FERS_Board{self.boardNo}_energyHG_{self.channelNo}"

    def GetLGChannelName(self):
        return f"FERS_Board{self.boardNo}_energyLG_{self.channelNo}"


class DRSChannel(CaloXChannel):
    def __init__(self, iTowerX: int, iTowerY: int, iBoardX: int, iBoardY: int, isCer: bool, channelNo: int, groupNo: int, boardNo: int):
        super().__init__(iTowerX, iTowerY, iBoardX, iBoardY, isCer)
        self.channelNo = channelNo
        self.groupNo = groupNo
        self.boardNo = boardNo

    def __str__(self):
        base_str = super().__str__()
        base_str += (f", channelNo={self.channelNo}, groupNo={self.groupNo}, "
                     f"boardNo={self.boardNo})")
        base_str = base_str.replace("CaloXChannel", "DRSChannel")
        return base_str

    def __copy__(self):
        """
        Create a shallow copy of the DRS channel.
        """
        return DRSChannel(
            self.iTowerX, self.iTowerY, self.iBoardX, self.iBoardY,
            self.isCer, self.channelNo, self.groupNo, self.boardNo
        )

    def GetChannelName(self):
        return f"DRS_Board{self.boardNo}_Group{self.groupNo}_Channel{self.channelNo}"


class Board(object):
    """
    A class to represent a generic board.
    """

    def __init__(self, boardNo):
        self.boardNo = boardNo  # Board number
        self.channels = []  # List of channels on the board
        # channels are organized as 2D list of CaloXChannel objects
        # by (iBoardX, iBoardY)

    def __str__(self):
        s = f"Board Number: {self.boardNo}\n"
        for row in self.channels:
            for channel in row:
                s += str(channel) + "\n"
        return s.strip()

    def __iter__(self):
        for row in self.channels:
            for channel in row:
                yield channel

    def __getitem__(self, idx):
        """
        Allow index access to the board's channels.
        Usage: board[x, y] returns the channel at (x, y).
        """
        try:
            if isinstance(idx, tuple) and len(idx) == 2:
                iBoardX, iBoardY = idx
                return self.channels[iBoardX][iBoardY]
            else:
                return self.channels[idx]
        except IndexError:
            if len(self.channels) == 0 or len(self.channels[0]) == 0:
                raise IndexError("Board has no channels defined.")
            else:
                raise IndexError(
                    f"Channel not found, Range is (0-{len(self.channels)-1}, 0-{len(self.channels[0])-1})")

    def GetChannelByTower(self, iTowerX, iTowerY, isCer=False):
        """
        Get a channel by its tower coordinates (iTowerX, iTowerY).
        Returns the first matching channel or None if not found.
        """
        for row in self.channels:
            for channel in row:
                if channel.iTowerX == iTowerX and channel.iTowerY == iTowerY and channel.isCer == isCer:
                    return channel
        print(
            f"\033[91m Warning: Channel not found for tower ({iTowerX}, {iTowerY}) on board {self.boardNo}.\033[0m")
        return None

    def GetCerChannels(self):
        """
        Get all CER channels on the board.
        """
        cerchannels = []
        for row in self.channels:
            for channel in row:
                if channel.isCer:
                    cerchannels.append(channel)
        return cerchannels

    def GetSciChannels(self):
        """
        Get all SCI channels on the board.
        """
        scichannels = []
        for row in self.channels:
            for channel in row:
                if not channel.isCer:
                    scichannels.append(channel)
        return scichannels

    def GetListOfTowers(self):
        """
        Get a list of unique towers on the board.
        Returns a list of tuples (iTowerX, iTowerY).
        """
        towers = set()
        for row in self.channels:
            for channel in row:
                towers.add((channel.iTowerX, channel.iTowerY))
        return list(towers)

    def MapCerToSci(self):
        """
        Map CER channels to SCI channels.
        Returns a dictionary mapping CER channels to SCI channels.
        """
        cer_channels = self.GetCerChannels()
        sci_channels = self.GetSciChannels()
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
        # organized by (iBoardX, iBoardY)
        self.channels = buildFERSBase(is6mm=is6mm, boardNo=boardNo)

    def __str__(self):
        base_str = super().__str__()
        base_str = "FERS " + base_str
        return base_str

    def __copy__(self):
        """
        Create a shallow copy of the FERS board.
        """
        new_board = FERSBoard(self.boardNo, self.is6mm)
        new_board.channels = [row.copy() for row in self.channels]
        return new_board

    def copy(self, boardNo):
        """
        Create a deep copy of the FERS board.
        """
        new_board = FERSBoard(boardNo, self.is6mm)
        new_board.channels = [[channel.__copy__() for channel in row]
                              for row in self.channels]
        for row in new_board.channels:
            for channel in row:
                channel.boardNo = boardNo
        return new_board

    def Is6mm(self):
        return self.is6mm

    def Is3mm(self):
        return not self.Is6mm()


class DRSBoard(Board):
    """
    A class to represent a DRS board.
    """

    def __init__(self, boardNo):
        """
        Initialize a DRS board with a specific board number and channel configuration.
        :param boardNo: The board number.
        """
        super().__init__(boardNo)
        # channels is a 2D list of DRSChannel objects
        # organized by (iBoardX, iBoardY)
        self.channels = buildDRSBase(boardNo=boardNo)

    def __str__(self):
        base_str = super().__str__()
        base_str = "DRS " + base_str
        return base_str

    def __copy__(self):
        """
        Create a shallow copy of the DRS board.
        """
        new_board = DRSBoard(self.boardNo)
        new_board.channels = [row.copy() for row in self.channels]
        return new_board

    def copy(self, boardNo):
        """
        Create a deep copy of the DRS board.
        """
        new_board = DRSBoard(boardNo)
        new_board.channels = [[channel.__copy__() for channel in row]
                              for row in self.channels]
        for row in new_board.channels:
            for channel in row:
                channel.boardNo = boardNo
        return new_board


# physical channels to the readout channel names
A5202_map_tmp = np.array([
    [0, 2, 1, 3],
    [4, 6, 5, 7],
    [8, 10, 9, 11],
    [12, 14, 13, 15],
    [16, 18, 17, 19],
    [20, 22, 21, 23],
    [24, 26, 25, 27],
    [28, 30, 29, 31],
    [32, 34, 33, 35],
    [36, 38, 37, 39],
    [40, 42, 41, 43],
    [44, 46, 45, 47],
    [48, 50, 49, 51],
    [52, 54, 53, 55],
    [56, 58, 57, 59],
    [60, 62, 61, 63]
], dtype=int)
# swap the axes such that iX and iY are correct
A5202_map = np.swapaxes(A5202_map_tmp, 0, 1)

drs_map_tmp = np.array([
    [3, 2, 1, 0],
    [7, 6, 5, 4],
    [11, 10, 9, 8],
    [15, 14, 13, 12],
    [19, 18, 17, 16],
    [23, 22, 21, 20],
    [27, 26, 25, 24],
    [31, 30, 29, 28],
], dtype=int)
drs_map = np.swapaxes(drs_map_tmp, 0, 1)


def buildFERSBase(is6mm=False, boardNo=0):
    channels_FERS = []
    for ix in range(0, 4):
        channels_FERS_one_row = []
        for iy in range(0, 16):
            channelNo = A5202_map[ix, iy]
            if is6mm:
                isCer = (iy % 2 == 0)
                channel = FERSChannel(ix, -int(iy/2), ix, iy, isCer,
                                      channelNo, boardNo)
                channels_FERS_one_row.append(channel)
            else:
                isCer = (ix % 2 == 1)
                # 3mm has higher granularity in global Y
                channel = FERSChannel(
                    int(ix/2), -float(iy)/4, ix, iy, isCer, channelNo, boardNo)
                channels_FERS_one_row.append(channel)
        channels_FERS.append(channels_FERS_one_row)
    return channels_FERS


def buildDRSBase(boardNo=0):
    # right now only have DRS on 6mm
    channels_DRS = []
    for ix in range(0, 4):
        channels_DRS_one_row = []
        for iy in range(0, 8):
            channelNo = drs_map[ix, iy]
            isCer = (iy % 2 == 0)
            # DRS channels are grouped in groups of 8
            # groupNo is the group number (0-3)
            # chanNo is the channel number within the group (0-7)
            groupNo = (channelNo // 8)
            chanNo = (channelNo % 8)
            channel = DRSChannel(ix, -int(iy/2), ix, iy, isCer,
                                 chanNo, groupNo, boardNo)
            channels_DRS_one_row.append(channel)
        channels_DRS.append(channels_DRS_one_row)
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

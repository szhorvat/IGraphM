
void message(const char *msg, MessageType type) {
    if (msg == NULL)
        return;

    if (libData->AbortQ())
        return; // trying to use the MathLink connection during an abort will break it

    const char *tag;
    switch (type) {
    case M_ERROR:
        tag = "error";
        break;
    case M_WARNING:
        tag = "warning";
        break;
    case M_ASSERT:
        tag = "assert";
        break;
    case M_INFO:
    default:
        tag = "info";
    }

    MLINK link = libData->getMathLink(libData);
    MLPutFunction(link, "EvaluatePacket", 1);
        MLPutFunction(link, "Message", 2);
            MLPutFunction(link, "MessageName", 2);
                MLPutSymbol(link, LTEMPLATE_MESSAGE_SYMBOL);
                MLPutString(link, tag);
            MLPutString(link, msg);
    libData->processMathLink(link);
    int pkt = MLNextPacket(link);
    if (pkt == RETURNPKT)
        MLNewPacket(link);
}

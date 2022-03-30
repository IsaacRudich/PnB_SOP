"""
    sequencingRRB(setting::RRBSetting, node::RestrictedSequencingNode, packet)

Call the correct rrb method based on the settings

# Arguments
- `setting::RRBSetting`: Which function to call
- `node::RestrictedSequencingNode`: The current node
- `packet`: Data with the necesary info
"""
function sequencingRRB(setting::RRBSetting, node::RestrictedSequencingNode, packet)
    if setting == sop
        return sopRRB(node.visited, getState(node), node.value, packet[1],packet[2])
    end
end
